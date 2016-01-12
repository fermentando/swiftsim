/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "task.h"

/* Local headers. */
#include "atomic.h"
#include "error.h"
#include "lock.h"

struct task *store;

/* Task type names. */
const char *taskID_names[task_type_count] = {
    "none",    "sort",    "self",    "pair",     "sub",  "init",
    "ghost",   "drift",   "kick",    "send",     "recv",
    "grav_pp", "grav_mm", "grav_up", "grav_down",
    "psort", "split_cell", "rewait"};

const char *subtaskID_names[task_type_count] = {"none",  "density",
                                                "force", "grav"};

/**
 * @brief Unlock the cell held by this task.
 *
 * @param t The #task.
 */

void task_unlock(struct task *t) {

  /* Act based on task type. */
  switch (t->type) {
    case task_type_self:
    case task_type_sort:
      cell_unlocktree(t->ci);
      break;
    case task_type_pair:
    case task_type_sub:
      cell_unlocktree(t->ci);
      if (t->cj != NULL) cell_unlocktree(t->cj);
      break;
    case task_type_grav_pp:
    case task_type_grav_mm:
    case task_type_grav_down:
      cell_gunlocktree(t->ci);
      if (t->cj != NULL) cell_gunlocktree(t->cj);
      break;
    default:
      break;
  }
}

/**
 * @brief Try to lock the cells associated with this task.
 *
 * @param t the #task.
 */

int task_lock(struct task *t) {

  int type = t->type;
  struct cell *ci = t->ci, *cj = t->cj;

  /* Communication task? */
  if (type == task_type_recv || type == task_type_send) {

#ifdef WITH_MPI
    /* Check the status of the MPI request. */
    int res, err;
    MPI_Status stat;
    if ((err = MPI_Test(&t->req, &res, &stat)) != MPI_SUCCESS) {
      char buff[MPI_MAX_ERROR_STRING];
      int len;
      MPI_Error_string(err, buff, &len);
      error("Failed to test request on send/recv task (tag=%i, %s).", t->flags,
            buff);
    }
    return res;
#else
    error("SWIFT was not compiled with MPI support.");
#endif

  }

  /* Unary lock? */
  else if (type == task_type_self || type == task_type_sort ||
           (type == task_type_sub && cj == NULL)) {
    if (cell_locktree(ci) != 0) return 0;
  }

  /* Otherwise, binary lock. */
  else if (type == task_type_pair || (type == task_type_sub && cj != NULL)) {
    if (ci->hold || cj->hold) return 0;
    if (cell_locktree(ci) != 0) return 0;
    if (cell_locktree(cj) != 0) {
      cell_unlocktree(ci);
      return 0;
    }
  }

  /* Gravity tasks? */
  else if (type == task_type_grav_mm || type == task_type_grav_pp ||
           type == task_type_grav_down) {
    if (ci->ghold || (cj != NULL && cj->ghold)) return 0;
    if (cell_glocktree(ci) != 0) return 0;
    if (cj != NULL && cell_glocktree(cj) != 0) {
      cell_gunlocktree(ci);
      return 0;
    }
  }

  /* If we made it this far, we've got a lock. */
  return 1;
}

/**
 * @brief Remove all unlocks to tasks that are of the given type.
 *
 * @param t The #task.
 * @param type The task type ID to remove.
 */

void task_cleanunlock(struct task *t, int type) {

  int k;

  lock_lock(&t->lock);

  for (k = 0; k < t->nr_unlock_tasks; k++)
    if (t->unlock_tasks[k]->type == type) {
      t->nr_unlock_tasks -= 1;
      t->unlock_tasks[k] = t->unlock_tasks[t->nr_unlock_tasks];
    }

  lock_unlock_blind(&t->lock);
}

/**
 * @brief Remove an unlock_task from the given task.
 *
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 */

void task_rmunlock(struct task *ta, struct task *tb) {

  int k;

  lock_lock(&ta->lock);

  for (k = 0; k < ta->nr_unlock_tasks; k++)
    if (ta->unlock_tasks[k] == tb) {
      ta->nr_unlock_tasks -= 1;
      ta->unlock_tasks[k] = ta->unlock_tasks[ta->nr_unlock_tasks];
      lock_unlock_blind(&ta->lock);
      return;
    }
  error("Task not found.");
}

/**
 * @brief Remove an unlock_task from the given task.
 *
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 *
 * Differs from #task_rmunlock in that it will not fail if
 * the task @c tb is not in the unlocks of @c ta.
 */

void task_rmunlock_blind(struct task *ta, struct task *tb) {

  int k;

  lock_lock(&ta->lock);

  for (k = 0; k < ta->nr_unlock_tasks; k++)
    if (ta->unlock_tasks[k] == tb) {
      ta->nr_unlock_tasks -= 1;
      ta->unlock_tasks[k] = ta->unlock_tasks[ta->nr_unlock_tasks];
      break;
    }

  lock_unlock_blind(&ta->lock);
}

/**
 * @brief Add an unlock_task to the given task.
 *
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 */

void task_addunlock(struct task *ta, struct task *tb) {

  error("Use sched_addunlock instead.");

  /* Add the lock atomically. */
  ta->unlock_tasks[atomic_inc(&ta->nr_unlock_tasks)] = tb;

  /* Check a posteriori if we did not overshoot. */
  if (ta->nr_unlock_tasks > task_maxunlock)
    error("Too many unlock_tasks in task.");
}

void task_addunlock_old(struct task *ta, struct task *tb) {

  int k;

  lock_lock(&ta->lock);

  /* Check if ta already unlocks tb. */
  for (k = 0; k < ta->nr_unlock_tasks; k++)
    if (ta->unlock_tasks[k] == tb) {
      error("Duplicate unlock.");
      lock_unlock_blind(&ta->lock);
      return;
    }

  if (ta->nr_unlock_tasks == task_maxunlock)
    error("Too many unlock_tasks in task.");

  ta->unlock_tasks[ta->nr_unlock_tasks] = tb;
  ta->nr_unlock_tasks += 1;

  lock_unlock_blind(&ta->lock);
}


void task_print_mask(unsigned int mask) {

  int k;

  printf("task_print_mask: The tasks to run are [");
  for (k = 1; k < task_type_count; k++)
    printf(" %s=%s", taskID_names[k], (mask & (1 << k)) ? "yes" : "no");
  printf(" ]\n");
}


void task_do_rewait(struct task *t) {

	  const unsigned int mask = t->flags;
	  
	  for (struct task *t2 = (struct task *)t->ci; t2 != (struct task *)t->cj; t2++) {

	    if( t2->skip ) continue;
	    
	    /* Skip tasks not in the mask */
	    if( !((1<<t2->type) & mask) ) continue;

	    /* Skip sort tasks that have already been performed */
	    if(t2->type == task_type_sort && t2->flags == 0) continue;

	    if(store == NULL && t2->type==task_type_pair && t2->subtype==task_subtype_density) {
	      message("\n");
	      message("Checking task %s-%s address: %p", taskID_names[t2->type], subtaskID_names[t2->subtype], t2);
	      store = t2;
	    }

	    for (int k = 0; k < t2->nr_unlock_tasks; k++) {

	      struct task *t3=t2->unlock_tasks[k];
	      
	      atomic_inc(&t3->wait);

	      if (t3 == store) {
		message("Unlocked by task %s-%s address: %p" , taskID_names[t2->type], subtaskID_names[t2->subtype], t2);
	      }
	      
	    }
	  }


  
}
