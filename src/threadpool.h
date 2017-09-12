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
#ifndef SWIFT_THREADPOOL_H
#define SWIFT_THREADPOOL_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <pthread.h>
#include <semaphore.h>

/* Local includes. */
#include "cycle.h"

/* Local defines. */
#define threadpool_log_initial_size 1000
#define threadpool_default_chunk_ratio 7

/* Forward declaration of the threadpool struct. */
struct threadpool;

/* Function type for mappings. */
typedef void (*threadpool_map_function)(void *map_data, int num_elements,
                                        void *extra_data);
typedef void (*threadpool_rmap_function)(struct threadpool *t, void *map_data,
                                         void *extra_data);

/* Data for threadpool logging. */
struct mapper_log_entry {

  /* ID of the thread executing the chunk. */
  int tid;

  /* Size of the chunk processed. */
  int chunk_size;

  /* Pointer to the mapper function. */
  void *map_function;

  /*! Start and end time of this task */
  ticks tic, toc;
};

struct mapper_log {
  /* Log of threadpool mapper calls. */
  struct mapper_log_entry *log;

  /* Size of the allocated log. */
  int size;

  /* Number of entries in the log. */
  int count;
};

/* Enum for threadpool mode of operation. */
enum threadpool_mode {
  threadpool_mode_none,
  threadpool_mode_map,
  threadpool_mode_rmap
};

/* Data of a threadpool. */
struct threadpool {

  /* The threads themselves. */
  pthread_t *threads;

  /* This is where threads go to rest. */
  pthread_barrier_t wait_barrier;
  pthread_barrier_t run_barrier;

  /* Current mode of operation. */
  enum threadpool_mode mode;

  /* Current map data and count. */
  void *map_data, *map_extra_data;
  volatile size_t map_data_count, map_data_size, map_data_stride,
      map_data_chunk;
  volatile threadpool_map_function map_function;

  /* Re-entrant mapping data. */
  volatile void *volatile *rmap_data;
  size_t rmap_data_size;
  void *rmap_extra_data;
  volatile size_t rmap_first, rmap_last;
  volatile size_t rmap_waiting;
  volatile threadpool_rmap_function rmap_function;
  sem_t rmap_semaphore;

  /* Number of threads in this pool. */
  int num_threads;

  /* Counter for the number of threads that are done. */
  volatile int num_threads_running;

#ifdef SWIFT_DEBUG_THREADPOOL
  struct mapper_log *logs;
#endif
};

/* Function prototypes. */
void threadpool_init(struct threadpool *tp, int num_threads);
void threadpool_map(struct threadpool *tp, threadpool_map_function map_function,
                    void *map_data, size_t N, int stride, int chunk,
                    void *extra_data);
void threadpool_rmap_add(struct threadpool *tp, void **map_data, size_t count);
void threadpool_rmap(struct threadpool *tp,
                     threadpool_rmap_function rmap_function, void **map_data,
                     size_t count, size_t size, void *extra_data);
void threadpool_clean(struct threadpool *tp);
#ifdef SWIFT_DEBUG_THREADPOOL
void threadpool_reset_log(struct threadpool *tp);
void threadpool_dump_log(struct threadpool *tp, const char *filename,
                         int reset);
#endif

#endif /* SWIFT_THREADPOOL_H */
