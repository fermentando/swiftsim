/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Config parameters. */
#include "../config.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "mpi_mesh_gravity.h"

/* Local includes. */
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "part.h"
#include "space.h"

/**
 * @brief Increment the value associated with a hashmap key
 *
 * The hashmap entry specified by key is incremented by m_add
 * if it already exists, otherwise it is created and set equal
 * to m_add.
 *
 * @param map Pointer to the hash map.
 * @param key Which hash key to update.
 * @param m_add Amount by which to increment the value in the hash map.
 */
__attribute__((always_inline, const)) INLINE static void add_to_hashmap(
    hashmap_t *map, hashmap_key_t key, double m_add) {

  int created = 0;
  hashmap_value_t *value = hashmap_get_new(map, key, &created);
  if (created) {
    /* Key was not present, so this is a new element */
    value->value_dbl = m_add;
  } else {
    /* Key was present, so add m_add to previous value */
    value->value_dbl += m_add;
  }
}

/**
 * @brief Accumulate local contributions to the density field
 *
 * Creates a hashmap with the contributions to the density
 * mesh from local particles. Here we require that hash_key_t
 * can store values up to at least N*N*N.
 *
 * TODO: parallelize. E.g. one hashmap per thread and combine
 * when done.
 *
 * @param N The size of the mesh
 * @param fac Inverse of the cell size
 * @param s The #space containing the particles.
 * @param map The hashmap in which to store the results
 *
 */
void accumulate_local_gparts_to_hashmap(const int N, const double fac,
                                        const struct space *s, hashmap_t *map) {

  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int *local_cells = s->local_cells_top;
  const int nr_local_cells = s->nr_local_cells;

  /* Loop over local cells */
  for (int icell = 0; icell < nr_local_cells; icell += 1) {
    /* Get a pointer to this cell */
    struct cell *cell = &(s->cells_top[local_cells[icell]]);
    /* Loop over particles in this cell */
    for (int ipart = 0; ipart < cell->grav.count; ipart += 1) {

      struct gpart *gp = &(cell->grav.parts[ipart]);

      /* Box wrap the multipole's position */
      const double pos_x = box_wrap(gp->x[0], 0., dim[0]);
      const double pos_y = box_wrap(gp->x[1], 0., dim[1]);
      const double pos_z = box_wrap(gp->x[2], 0., dim[2]);

      /* Workout the CIC coefficients */
      int i = (int)(fac * pos_x);
      if (i >= N) i = N - 1;
      const double dx = fac * pos_x - i;
      const double tx = 1. - dx;

      int j = (int)(fac * pos_y);
      if (j >= N) j = N - 1;
      const double dy = fac * pos_y - j;
      const double ty = 1. - dy;

      int k = (int)(fac * pos_z);
      if (k >= N) k = N - 1;
      const double dz = fac * pos_z - k;
      const double tz = 1. - dz;

#ifdef SWIFT_DEBUG_CHECKS
      if (i < 0 || i >= N) error("Invalid gpart position in x");
      if (j < 0 || j >= N) error("Invalid gpart position in y");
      if (k < 0 || k >= N) error("Invalid gpart position in z");
#endif

      /* Accumulate contributions to the hashmap */
      const double mass = gp->mass;
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t(i + 0, j + 0, k + 0, N),
          mass * tx * ty * tz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t(i + 0, j + 0, k + 1, N),
          mass * tx * ty * dz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t(i + 0, j + 1, k + 0, N),
          mass * tx * dy * tz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t(i + 0, j + 1, k + 1, N),
          mass * tx * dy * dz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t(i + 1, j + 0, k + 0, N),
          mass * dx * ty * tz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t(i + 1, j + 0, k + 1, N),
          mass * dx * ty * dz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t(i + 1, j + 1, k + 0, N),
          mass * dx * dy * tz);
      add_to_hashmap(
          map,
          (hashmap_key_t)row_major_id_periodic_size_t(i + 1, j + 1, k + 1, N),
          mass * dx * dy * dz);

    } /* Next particle */
  }   /* Next cell */

  return;
}

/**
 * @brief Store contributions to the mesh as (index, mass) pairs
 */
struct mesh_key_value {
  hashmap_key_t key;
  double value;
};

/**
 * @brief Comparison function to sort mesh_key_value by key
 *
 * @param a The first #mesh_key_value object.
 * @param b The second #mesh_key_value object.
 * @return 1 if a's key field is greater than b's key field,
 * -1 if a's key field is less than b's key field, and zero otherwise
 *
 */
int cmp_func_mesh_key_value(const void *a, const void *b) {
  struct mesh_key_value *a_key_value = (struct mesh_key_value *)a;
  struct mesh_key_value *b_key_value = (struct mesh_key_value *)b;
  if (a_key_value->key > b_key_value->key)
    return 1;
  else if (a_key_value->key < b_key_value->key)
    return -1;
  else
    return 0;
}

/**
 * @brief Data needed to iterate over a hashmap copying key-value pairs
 */
struct hashmap_mapper_data {
  size_t offset;
  struct mesh_key_value *buf;
};

/**
 * @brief Mapper function to copy elements from the hashmap
 *
 * @param key The key associated with this hashmap entry
 * @param value The value associated with this hashmap entry
 * @param data Contains pointer to output buffer and offset to
 * next element to write
 *
 */
void hashmap_copy_elements_mapper(hashmap_key_t key, hashmap_value_t *value,
                                  void *data) {

  struct hashmap_mapper_data *mapper_data = (struct hashmap_mapper_data *)data;
  struct mesh_key_value *element = &(mapper_data->buf[mapper_data->offset]);
  element->key = key;
  element->value = value->value_dbl;
  mapper_data->offset += 1;
}

/**
 * @brief Given an array of struct mesh_key_value, send nr_send[i]
 * elements to each node i. Allocates the receive buffer recvbuf
 * to the appropriate size and returns its size in nr_recv_tot.
 *
 * @param nr_send Number of elements to send to each other node
 * @param sendbuf The elements to send
 * @param nr_recv_tot Returns total number of elements received
 * @param recvbuf Returns a pointer to the newly received data
 *
 */
void exchange_mesh_cells(size_t *nr_send, struct mesh_key_value *sendbuf,
                         size_t *nr_recv_tot, struct mesh_key_value **recvbuf) {

  /* Determine rank, number of ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  /* Determine number of elements to receive */
  size_t *nr_recv = malloc(nr_nodes * sizeof(size_t));
  MPI_Alltoall(nr_send, sizeof(size_t), MPI_BYTE, nr_recv, sizeof(size_t),
               MPI_BYTE, MPI_COMM_WORLD);

  /* Find total elements to receive */
  *nr_recv_tot = 0;
  for (int i = 0; i < nr_nodes; i += 1) {
    *nr_recv_tot += nr_recv[i];
  }

  /* Allocate the receive buffer */
  if (swift_memalign("mesh_recvbuf", (void **)recvbuf, 32,
                     *nr_recv_tot * sizeof(struct mesh_key_value)) != 0)
    error("Failed to allocate receive buffer for constructing MPI FFT mesh");

  /* Compute send offsets */
  size_t *send_offset = malloc(nr_nodes * sizeof(size_t));
  send_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    send_offset[i] = send_offset[i - 1] + nr_send[i - 1];
  }

  /* Compute receive offsets */
  size_t *recv_offset = malloc(nr_nodes * sizeof(size_t));
  recv_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    recv_offset[i] = recv_offset[i - 1] + nr_recv[i - 1];
  }

  /* Allocate request objects (one send and receive per node) */
  MPI_Request *request = malloc(2 * sizeof(MPI_Request) * nr_nodes);

  /* Make type to communicate mesh_key_value struct */
  MPI_Datatype mesh_key_value_mpi_type;
  if (MPI_Type_contiguous(sizeof(struct mesh_key_value) / sizeof(unsigned char),
                          MPI_BYTE, &mesh_key_value_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&mesh_key_value_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for mesh_key_value struct.");
  }

  /*
   * Post the send operations. This is an alltoallv really but
   * we want to avoid the limits imposed by int counts and offsets
   * in MPI_Alltoallv.
   */
  for (int i = 0; i < nr_nodes; i += 1) {
    if (nr_send[i] > 0) {
      MPI_Isend(&(sendbuf[send_offset[i]]), (int)nr_send[i],
                mesh_key_value_mpi_type, i, 0, MPI_COMM_WORLD, &(request[i]));
    } else {
      request[i] = MPI_REQUEST_NULL;
    }
  }

  /* Post the receives */
  for (int i = 0; i < nr_nodes; i += 1) {
    if (nr_recv[i] > 0) {
      MPI_Irecv(&((*recvbuf)[recv_offset[i]]), (int)nr_recv[i],
                mesh_key_value_mpi_type, i, 0, MPI_COMM_WORLD,
                &(request[i + nr_nodes]));
    } else {
      request[i + nr_nodes] = MPI_REQUEST_NULL;
    }
  }

  /* Wait for everything to complete */
  MPI_Waitall(2 * nr_nodes, request, MPI_STATUSES_IGNORE);

  /* Done with the MPI type */
  MPI_Type_free(&mesh_key_value_mpi_type);

  /* Tidy up */
  free(recv_offset);
  free(send_offset);
  free(request);
  free(nr_recv);
}

/**
 * @brief Convert hashmaps to a slab-distributed 3D mesh
 *
 * For FFTW each rank needs to hold a slice of the full mesh.
 * This routine does the necessary communication to convert
 * the per-rank hashmaps into a slab-distributed mesh.
 *
 * @param e Pointer to the engine struct
 * @param N The size of the mesh
 * @param Nslice The thickness of the slice to store on this rank
 * @param map The hashmap with the local part of the mesh
 * @param mesh Pointer to the output data buffer
 *
 */
void hashmaps_to_slices(const int N, const int Nslice, hashmap_t *map,
                        double *mesh) {

  /* Determine rank, number of ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  /* Allocate storage for mesh cells to send */
  size_t nr_send_tot = hashmap_size(map);
  struct mesh_key_value *mesh_sendbuf;
  if (swift_memalign("mesh_sendbuf", (void **)&mesh_sendbuf, 32,
                     nr_send_tot * sizeof(struct mesh_key_value)) != 0)
    error("Failed to allocate send buffer for constructing MPI FFT mesh");

  /* Copy the key-value pairs to the new array */
  struct hashmap_mapper_data mapper_data = {(size_t)0, mesh_sendbuf};
  hashmap_iterate(map, hashmap_copy_elements_mapper, &mapper_data);

  /* Then sort the array of local mesh cells by key. This means
   * they're sorted by x coordinate, then y coordinate, then z coordinate.
   * We're going to distribute them between ranks according to their
   * x coordinate, so this puts them in order of destination rank.
   */
  qsort(mesh_sendbuf, nr_send_tot, sizeof(struct mesh_key_value),
        cmp_func_mesh_key_value);

  /* Get width of the slice on each rank */
  int *slice_width = malloc(sizeof(int) * nr_nodes);
  MPI_Allgather(&Nslice, 1, MPI_INT, slice_width, 1, MPI_INT, MPI_COMM_WORLD);

  /* Determine offset to the slice on each rank */
  int *slice_offset = malloc(sizeof(int) * nr_nodes);
  slice_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    slice_offset[i] = slice_offset[i - 1] + slice_width[i - 1];
  }

  /* Compute how many elements are to be sent to each rank */
  size_t *nr_send = malloc(nr_nodes * sizeof(size_t));
  for (int i = 0; i < nr_nodes; i += 1) {
    nr_send[i] = 0;
  }
  int dest_node = 0;
  for (size_t i = 0; i < nr_send_tot; i += 1) {
    /* Get the x coordinate of this mesh cell */
    int mesh_x = mesh_sendbuf[i].key / ((hashmap_key_t)(N * N));
    /* Advance to the destination node that is to contain this x coordinate */
    while ((mesh_x >= slice_offset[dest_node] + slice_width[dest_node]) ||
           (slice_width[dest_node] == 0)) {
      dest_node += 1;
    }
    nr_send[dest_node] += 1;
  }

  /* Carry out the communication */
  size_t nr_recv_tot;
  struct mesh_key_value *mesh_recvbuf;
  exchange_mesh_cells(nr_send, mesh_sendbuf, &nr_recv_tot, &mesh_recvbuf);

  /* No longer need the send buffer */
  swift_free("mesh_sendbuf", mesh_sendbuf);
  free(nr_send);

  /* Copy received data to the output buffer */
  for (size_t i = 0; i < nr_recv_tot; i += 1) {
#ifdef SWIFT_DEBUG_CHECKS
    if (mesh_recvbuf[i].key < N * N * slice_offset[nodeID])
      error("Received cell's local index is negative");
    if (mesh_recvbuf[i].key - N * N * slice_offset[nodeID] >=
        N * N * slice_width[nodeID])
      error("Received cell's local index is too large");
#endif
    mesh[mesh_recvbuf[i].key - (N * N * slice_offset[nodeID])] +=
        mesh_recvbuf[i].value;
  }

  /* Tidy up */
  free(slice_width);
  free(slice_offset);
  swift_free("mesh_recvbuf", mesh_recvbuf);
}
