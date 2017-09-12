/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef SWIFT_DEBUG_THREADPOOL
#include <dlfcn.h>
#endif

/* This object's header. */
#include "threadpool.h"

/* Local headers. */
#include "atomic.h"
#include "clocks.h"
#include "error.h"
#include "minmax.h"

#ifdef SWIFT_DEBUG_THREADPOOL
/**
 * @brief Store a log entry of the given chunk.
 */
void threadpool_log(struct threadpool *tp, int tid, size_t chunk_size,
                    ticks tic, ticks toc) {
  struct mapper_log *log = &tp->logs[tid > 0 ? tid : 0];

  /* Check if we need to re-allocate the log buffer. */
  if (log->count == log->size) {
    log->size *= 2;
    struct mapper_log_entry *new_log;
    if ((new_log = (struct mapper_log_entry *)malloc(
             sizeof(struct mapper_log_entry) * log->size)) == NULL)
      error("Failed to re-allocate mapper log.");
    memcpy(new_log, log->log, sizeof(struct mapper_log_entry) * log->count);
    free(log->log);
    log->log = new_log;
  }

  /* Store the new entry. */
  struct mapper_log_entry *entry = &log->log[log->count];
  entry->tid = tid;
  entry->chunk_size = chunk_size;
  entry->tic = tic;
  entry->toc = toc;
  switch (tp->mode) {
    case threadpool_mode_map:
      entry->map_function = tp->map_function;
      break;
    case threadpool_mode_rmap:
      entry->map_function = tp->rmap_function;
      break;
    default:
      error("Can't log when not running.");
  }
  log->count++;
}

void threadpool_dump_log(struct threadpool *tp, const char *filename,
                         int reset) {

  /* Open the output file. */
  FILE *fd;
  if ((fd = fopen(filename, "w")) == NULL)
    error("Failed to create log file '%s'.", filename);

  /* Create a buffer of function names. */
  const int max_names = 100;
  struct name_entry {
    threadpool_map_function map_function;
    const char *name;
  };
  struct name_entry names[max_names];
  bzero(names, sizeof(struct name_entry) * max_names);

  /* Write a header. */
  fprintf(fd, "# map_function thread_id chunk_size tic toc\n");
  fprintf(fd, "# {'num_threads': %i, 'cpufreq': %lli}\n", tp->num_threads,
          clocks_get_cpufreq());

  /* Loop over the per-tid logs and dump them. */
  for (int k = 0; k < tp->num_threads; k++) {
    struct mapper_log *log = &tp->logs[k];

    /* Loop over the log entries and dump them. */
    for (int i = 0; i < log->count; i++) {

      struct mapper_log_entry *entry = &log->log[i];

      /* Look for the function pointer in the buffer. */
      int nid = 0;
      while (nid < max_names && names[nid].map_function != entry->map_function)
        nid++;

      /* If the name was not found, make a new entry. */
      if (nid == max_names) {
        for (int j = 1; j < max_names; j++) names[j - 1] = names[j];
        names[0].map_function = entry->map_function;
        Dl_info dl_info;
        dladdr(entry->map_function, &dl_info);
        names[0].name = dl_info.dli_sname;
        nid = 0;
      }

      /* Log a line to the file. */
      fprintf(fd, "%s %i %i %lli %lli\n", names[nid].name, entry->tid,
              entry->chunk_size, entry->tic, entry->toc);
    }

    /* Clear the log if requested. */
    if (reset) log->count = 0;
  }

  /* Close the file. */
  fclose(fd);
}
#endif  // SWIFT_DEBUG_THREADPOOL

/**
 * @brief Runner main loop, get a chunk and call the mapper function.
 */
void threadpool_chomp(struct threadpool *tp, int tid) {

  /* Loop until we can't get a chunk. */
  while (1) {
    /* Desired chunk size. */
    size_t chunk_size =
        (tp->map_data_size - tp->map_data_count) / (2 * tp->num_threads);
    if (chunk_size > tp->map_data_chunk) chunk_size = tp->map_data_chunk;
    if (chunk_size < 1) chunk_size = 1;

    /* Get a chunk and check its size. */
    size_t task_ind = atomic_add(&tp->map_data_count, chunk_size);
    if (task_ind >= tp->map_data_size) break;
    if (task_ind + chunk_size > tp->map_data_size)
      chunk_size = tp->map_data_size - task_ind;

/* Call the mapper function. */
#ifdef SWIFT_DEBUG_THREADPOOL
    ticks tic = getticks();
#endif
    tp->map_function((char *)tp->map_data + (tp->map_data_stride * task_ind),
                     chunk_size, tp->map_extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
    threadpool_log(tp, tid, chunk_size, tic, getticks());
#endif
  }
}

/**
 * @brief Runner main re-entrant loop.
 */
void threadpool_rchomp(struct threadpool *tp, int tid) {

  /* Main loop. */
  while (1) {

    /* Wait at the semaphore. */
    sem_wait(&tp->rmap_semaphore);
    if (tp->rmap_waiting == 0) return;

    /* Try to get an entry. */
    size_t ind = atomic_inc(&tp->rmap_first) % tp->rmap_data_size;
    while (tp->rmap_data[ind] == NULL) {
      if (tp->rmap_waiting == 0) return;
    }

    /* Yay, we got an element! */
    void *data = (void *)tp->rmap_data[ind];
    tp->rmap_data[ind] = NULL;

/* Call the mapper function. */
#ifdef SWIFT_DEBUG_THREADPOOL
    ticks tic = getticks();
#endif
    tp->rmap_function(tp, data, tp->rmap_extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
    threadpool_log(tp, tid, 1, tic, getticks());
#endif

    /* One less task waiting. If we were the last task, hit the
       semaphore to release the others. */
    if (atomic_dec(&tp->rmap_waiting) == 1) {
      for (int k = 0; k < tp->num_threads; k++) sem_post(&tp->rmap_semaphore);
      break;
    }
  }
}

void *threadpool_runner(void *data) {

  /* Our threadpool. */
  struct threadpool *tp = (struct threadpool *)data;

  /* Main loop. */
  while (1) {

    /* Let the controller know that this thread is waiting. */
    pthread_barrier_wait(&tp->wait_barrier);

    /* Wait for the controller. */
    pthread_barrier_wait(&tp->run_barrier);

    /* If no map function is specified, just die. We use this as a mechanism
       to shut down threads without leaving the barriers in an invalid state. */
    if (tp->map_function == NULL) pthread_exit(NULL);

    /* Do actual work. */
    switch (tp->mode) {
      case threadpool_mode_map:
        threadpool_chomp(tp, atomic_inc(&tp->num_threads_running));
        break;
      case threadpool_mode_rmap:
        threadpool_rchomp(tp, atomic_inc(&tp->num_threads_running));
        break;
      default:
        error("No mapping specified.");
    }
  }
}

/**
 * @brief Initialises the #threadpool with a given number of threads.
 *
 * @param tp The #threadpool.
 * @param num_threads The number of threads.
 */
void threadpool_init(struct threadpool *tp, int num_threads) {

  /* Initialize the thread counters. */
  tp->num_threads = num_threads;

#ifdef SWIFT_DEBUG_THREADPOOL
  if ((tp->logs = (struct mapper_log *)malloc(sizeof(struct mapper_log) *
                                              num_threads)) == NULL)
    error("Failed to allocate mapper logs.");
  for (int k = 0; k < num_threads; k++) {
    tp->logs[k].size = threadpool_log_initial_size;
    tp->logs[k].count = 0;
    if ((tp->logs[k].log = (struct mapper_log_entry *)malloc(
             sizeof(struct mapper_log_entry) * tp->logs[k].size)) == NULL)
      error("Failed to allocate mapper log.");
  }
#endif

  /* If there is only a single thread, do nothing more as of here as
     we will just do work in the (blocked) calling thread. */
  if (num_threads == 1) return;

  /* Init the barriers. */
  if (pthread_barrier_init(&tp->wait_barrier, NULL, num_threads) != 0 ||
      pthread_barrier_init(&tp->run_barrier, NULL, num_threads) != 0)
    error("Failed to initialize barriers.");

  /* Set the task counter to zero. */
  tp->map_data_size = 0;
  tp->map_data_count = 0;
  tp->map_data_stride = 0;
  tp->map_data_chunk = 0;
  tp->map_function = NULL;
  tp->rmap_function = NULL;
  tp->mode = threadpool_mode_none;

  /* Allocate the threads, one less than requested since the calling thread
     works as well. */
  if ((tp->threads = (pthread_t *)malloc(sizeof(pthread_t) *
                                         (num_threads - 1))) == NULL) {
    error("Failed to allocate thread array.");
  }

  /* Create and start the threads. */
  for (int k = 0; k < num_threads - 1; k++) {
    if (pthread_create(&tp->threads[k], NULL, &threadpool_runner, tp) != 0)
      error("Failed to create threadpool runner thread.");
  }

  /* Wait for all the threads to be up and running. */
  pthread_barrier_wait(&tp->wait_barrier);
}

/**
 * @brief Add a set of data pointers to a running re-entrant threadpool.
 *
 * @param tp The #threadpool, must be running as a re-entrant mapper.
 * @param map_data Pointer to an array of pointers containing the data to add.
 * @param count Number of elements in map_data.
 */
void threadpool_rmap_add(struct threadpool *tp, void **map_data, size_t count) {
  if (tp->mode != threadpool_mode_rmap)
    error("Threadpool not running as a re-entrant mapper.");

  /* Loop over the elements and add them to the rmap_data. */
  for (size_t k = 0; k < count; k++) {

    /* Skip empty elements. */
    if (map_data[k] == NULL) continue;

    /* Get an index in which to store the data. */
    size_t ind = atomic_inc(&tp->rmap_last) % tp->rmap_data_size;

    /* Store the data once the entry is clear. */
    atomic_inc(&tp->rmap_waiting);
    while (atomic_cas(&tp->rmap_data[ind], NULL, map_data[k]) != NULL)
      ;

    /* Post the item to the semaphore. */
    sem_post(&tp->rmap_semaphore);
  }
}

/**
 * @brief Map a function to an array of pointers where the function can add
 * elements to the list.
 *
 * @param tp The #threadpool on which to run.
 * @param rmap_function The function that will be applied to the map data.
 * @param map_data Array of pointers to the data, should be large enough to
 *        accommodate all in-flight re-entrant calls.
 * @param count Number of elements in map_data that have been set.
 * @param size Number of elements allocated in map_data.
 * @param extra_data Addtitional pointer that will be passed to the mapping
 *        function, may contain additional data.
 */
void threadpool_rmap(struct threadpool *tp,
                     threadpool_rmap_function rmap_function, void **map_data,
                     size_t count, size_t size, void *extra_data) {

#ifdef SWIFT_DEBUG_THREADPOOL
  ticks tic = getticks();
#endif

  /* Set up the rmap data. */
  tp->rmap_data = (volatile void **)map_data;
  tp->rmap_data_size = size;
  tp->rmap_first = 0;
  tp->rmap_last = count;
  tp->rmap_waiting = count;
  tp->rmap_function = rmap_function;
  tp->rmap_extra_data = extra_data;
  tp->num_threads_running = 0;
  tp->mode = threadpool_mode_rmap;

  /* Init the rmap semaphore with the initial number of items. */
  if (sem_init(&tp->rmap_semaphore, 0, count) != 0)
    error("Failed to initialize semaphore.");

  /* If it's just a single thread, don't bother with the barriers. */
  if (tp->num_threads == 1) {
    threadpool_rchomp(tp, 0);
  }

  /* Otherwise, launch all threads. */
  else {

    /* Wait for all the threads to be up and running. */
    pthread_barrier_wait(&tp->run_barrier);

    /* Do some work while I'm at it. */
    threadpool_rchomp(tp, tp->num_threads - 1);

    /* Wait for all threads to be done. */
    pthread_barrier_wait(&tp->wait_barrier);
  }
#ifdef SWIFT_DEBUG_THREADPOOL
  /* Log the total call time to thread id -1. */
  threadpool_log(tp, -1, tp->rmap_last, tic, getticks());
#endif

  /* Re-set the mode. */
  tp->mode = threadpool_mode_none;
}

/**
 * @brief Map a function to an array of data in parallel using a #threadpool.
 *
 * The function @c map_function is called on each element of @c map_data
 * in parallel.
 *
 * @param tp The #threadpool on which to run.
 * @param map_function The function that will be applied to the map data.
 * @param map_data The data on which the mapping function will be called.
 * @param N Number of elements in @c map_data.
 * @param stride Size, in bytes, of each element of @c map_data.
 * @param chunk Number of map data elements to pass to the function at a time,
 *        or zero to choose the number automatically.
 * @param extra_data Addtitional pointer that will be passed to the mapping
 *        function, may contain additional data.
 */
void threadpool_map(struct threadpool *tp, threadpool_map_function map_function,
                    void *map_data, size_t N, int stride, int chunk,
                    void *extra_data) {

#ifdef SWIFT_DEBUG_THREADPOOL
  ticks tic = getticks();
#endif

  /* If we just have a single thread, call the map function directly. */
  if (tp->num_threads == 1) {
    tp->mode = threadpool_mode_map;
    map_function(map_data, N, extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
    tp->map_function = map_function;
    threadpool_log(tp, 0, N, tic, getticks());
    threadpool_log(tp, -1, N, tic, getticks());
#endif
    tp->mode = threadpool_mode_none;
    return;
  }

  /* Set the map data and signal the threads. */
  tp->map_data_stride = stride;
  tp->map_data_size = N;
  tp->map_data_count = 0;
  tp->map_data_chunk =
      chunk ? chunk
            : max((int)(N / (tp->num_threads * threadpool_default_chunk_ratio)),
                  1);
  tp->map_function = map_function;
  tp->map_data = map_data;
  tp->map_extra_data = extra_data;
  tp->num_threads_running = 0;
  tp->mode = threadpool_mode_map;

  /* Wait for all the threads to be up and running. */
  pthread_barrier_wait(&tp->run_barrier);

  /* Do some work while I'm at it. */
  threadpool_chomp(tp, tp->num_threads - 1);

  /* Wait for all threads to be done. */
  pthread_barrier_wait(&tp->wait_barrier);

#ifdef SWIFT_DEBUG_THREADPOOL
  /* Log the total call time to thread id -1. */
  threadpool_log(tp, -1, N, tic, getticks());
#endif

  /* Re-set the mode. */
  tp->mode = threadpool_mode_none;
}

/**
 * @brief Re-sets the log for this #threadpool.
 */
#ifdef SWIFT_DEBUG_THREADPOOL
void threadpool_reset_log(struct threadpool *tp) {
  for (int k = 0; k < tp->num_threads; k++) tp->logs[k].count = 0;
}
#endif

/**
 * @brief Frees up the memory allocated for this #threadpool.
 */
void threadpool_clean(struct threadpool *tp) {

  if (tp->num_threads > 1) {
    /* Destroy the runner threads by calling them with a NULL mapper function
     * and waiting for all the threads to terminate. This ensures that no
     * thread is still waiting at a barrier. */
    tp->map_function = NULL;
    pthread_barrier_wait(&tp->run_barrier);
    for (int k = 0; k < tp->num_threads - 1; k++) {
      void *retval;
      pthread_join(tp->threads[k], &retval);
    }

    /* Release the barriers. */
    if (pthread_barrier_destroy(&tp->wait_barrier) != 0 ||
        pthread_barrier_destroy(&tp->run_barrier) != 0)
      error("Failed to destroy threadpool barriers.");

    /* Clean up memory. */
    free(tp->threads);
  }

#ifdef SWIFT_DEBUG_THREADPOOL
  for (int k = 0; k < tp->num_threads; k++) {
    free(tp->logs[k].log);
  }
  free(tp->logs);
#endif
}
