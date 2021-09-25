/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_DEBUGGING_GEAR_H
#define SWIFT_RT_DEBUGGING_GEAR_H

#ifdef SWIFT_RT_DEBUG_CHECKS

#include "rt_properties.h"

/**
 * @file src/rt/GEAR/rt_debugging.h
 * @brief Main header file for the GEAR radiative transfer scheme
 * extra debugging functions.
 */

/**
 * @brief Debugging checks loop over all star particles after each time step
 */
static void rt_debugging_end_of_step_stars_mapper(void *restrict map_data,
                                                  int scount,
                                                  void *restrict extra_data) {

  struct spart *restrict sparts = (struct spart *)map_data;
  const struct engine *restrict e = (struct engine *)extra_data;

  int emission_sum = 0;
  int iacts_with_parts_sum = 0;
  unsigned long long emission_sum_tot = 0ULL;
  unsigned long long iacts_with_parts_sum_tot = 0ULL;
  float emitted_energy[RT_NGROUPS];
  for (int g = 0; g < RT_NGROUPS; g++) emitted_energy[g] = 0.f;

  for (int k = 0; k < scount; k++) {

    struct spart *restrict sp = &sparts[k];
    emission_sum += sp->rt_data.debug_iact_hydro_inject;
    emission_sum_tot += sp->rt_data.debug_radiation_emitted_tot;
    /* Reset all values here in case stars won't be active next step */
    sp->rt_data.debug_iact_hydro_inject = 0;

    iacts_with_parts_sum += sp->rt_data.debug_iact_hydro_inject_prep;
    iacts_with_parts_sum_tot += sp->rt_data.debug_iact_hydro_inject_prep_tot;
    /* Reset all values here in case stars won't be active next step */
    sp->rt_data.debug_iact_hydro_inject_prep = 0;

    for (int g = 0; g < RT_NGROUPS; g++) {
      /* also check now that we actually injected the correct
       * amount of energy
       * sp->rt_data.emission_this_step: energy we should distribute
       *                                 this step
       * sp->rt_data.debug_injected_energy: energy we actually did
       *                                    distribute this step    */
      if (sp->rt_data.debug_injected_energy[g] != 0.f) {
        float diff = 1.f - sp->rt_data.emission_this_step[g] /
                               sp->rt_data.debug_injected_energy[g];

        if (fabs(diff) > 1e-3) {
          message(
              "Incorrect injection ID %lld: "
              "group %d expected %.3g got %.3g diff %.3g",
              sp->id, g, sp->rt_data.emission_this_step[g],
              sp->rt_data.debug_injected_energy[g], diff);
        }
      }
      emitted_energy[g] += sp->rt_data.debug_injected_energy[g];
    }

    for (int g = 0; g < RT_NGROUPS; g++) {
      sp->rt_data.debug_injected_energy[g] = 0.f;
    }
  }

  atomic_add(&e->rt_props->debug_radiation_emitted_this_step, emission_sum);
  atomic_add(&e->rt_props->debug_radiation_emitted_tot, emission_sum_tot);
  atomic_add(&e->rt_props->debug_star_injection_prep_iacts_with_parts_this_step, iacts_with_parts_sum);
  atomic_add(&e->rt_props->debug_star_injection_prep_iacts_with_parts_tot, iacts_with_parts_sum_tot);
  for (int g = 0; g < RT_NGROUPS; g++)
    atomic_add_f(&e->rt_props->debug_total_star_emitted_energy[g],
                 emitted_energy[g]);
}

/**
 * @brief Debugging checks loop over all hydro particles after each time step
 */
static void rt_debugging_end_of_step_hydro_mapper(void *restrict map_data,
                                                  int count,
                                                  void *restrict extra_data) {

  struct part *restrict parts = (struct part *)map_data;
  const struct engine *restrict e = (struct engine *)extra_data;

  int absorption_sum = 0;
  int iacts_with_stars_sum = 0;
  unsigned long long absorption_sum_tot = 0ULL;
  unsigned long long iacts_with_stars_sum_tot = 0ULL;
  float energy_sum[RT_NGROUPS];
  for (int g = 0; g < RT_NGROUPS; g++) energy_sum[g] = 0.f;
  float energy_density_sum[RT_NGROUPS];
  for (int g = 0; g < RT_NGROUPS; g++) energy_density_sum[g] = 0.f;

  for (int k = 0; k < count; k++) {

    struct part *restrict p = &parts[k];
    absorption_sum += p->rt_data.debug_iact_stars_inject;
    absorption_sum_tot += p->rt_data.debug_radiation_absorbed_tot;
    /* Reset all values here in case particles won't be active next step */
    p->rt_data.debug_iact_stars_inject = 0;

    iacts_with_stars_sum += p->rt_data.debug_iact_stars_inject_prep;
    iacts_with_stars_sum_tot += p->rt_data.debug_iact_stars_inject_prep_tot;
    /* Reset all values here in case particles won't be active next step */
    p->rt_data.debug_iact_stars_inject_prep = 0;

    for (int g = 0; g < RT_NGROUPS; g++) {
      energy_density_sum[g] += p->rt_data.density[g].energy;
      energy_sum[g] += p->rt_data.conserved[g].energy;
    }
  }

  atomic_add(&e->rt_props->debug_radiation_absorbed_this_step, absorption_sum);
  atomic_add(&e->rt_props->debug_radiation_absorbed_tot, absorption_sum_tot);
  atomic_add(&e->rt_props->debug_part_injection_prep_iacts_with_stars_this_step, iacts_with_stars_sum);
  atomic_add(&e->rt_props->debug_part_injection_prep_iacts_with_stars_tot, iacts_with_stars_sum_tot);

  for (int g = 0; g < RT_NGROUPS; g++) {
    atomic_add_f(&(e->rt_props->debug_total_radiation_conserved_energy[g]),
                 energy_sum[g]);
    atomic_add_f(&(e->rt_props->debug_total_radiation_energy_density[g]),
                 energy_density_sum[g]);
  }
}

/**
 * @brief At the end of each time step, loop over both hydro and star
 * particles and do whatever checks for this particular time step you
 * want done.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_checks_end_of_step(struct engine *e, int verbose) {

  struct space *s = e->s;
  if (!(e->policy & engine_policy_rt)) return;

  const ticks tic = getticks();

  /* reset values before the particle loops */
  e->rt_props->debug_radiation_emitted_this_step = 0;
  e->rt_props->debug_radiation_absorbed_this_step = 0;
  e->rt_props->debug_star_injection_prep_iacts_with_parts_this_step = 0;
  e->rt_props->debug_part_injection_prep_iacts_with_stars_this_step = 0;
  /* reset total counts as well. We track the totals since the beginning
   * of time in particles individually. */
  e->rt_props->debug_radiation_emitted_tot = 0LL;
  e->rt_props->debug_radiation_absorbed_tot = 0LL;
  e->rt_props->debug_star_injection_prep_iacts_with_parts_tot = 0LL;
  e->rt_props->debug_part_injection_prep_iacts_with_stars_tot = 0LL;
  for (int g = 0; g < RT_NGROUPS; g++) {
    e->rt_props->debug_total_radiation_energy_density[g] = 0.f;
    e->rt_props->debug_total_radiation_conserved_energy[g] = 0.f;
  }

  /* hydro particle loop */
  if (s->nr_parts > 0)
    threadpool_map(&e->threadpool, rt_debugging_end_of_step_hydro_mapper,
                   s->parts, s->nr_parts, sizeof(struct part),
                   threadpool_auto_chunk_size, /*extra_data=*/e);

  /* star particle loop */
  if (s->nr_sparts > 0)
    threadpool_map(&e->threadpool, rt_debugging_end_of_step_stars_mapper,
                   s->sparts, s->nr_sparts, sizeof(struct spart),
                   threadpool_auto_chunk_size, /*extra_data=*/e);

  /* message("This step:     %12d %12d %12d %12d", */
  /*           e->rt_props->debug_radiation_emitted_this_step, */
  /*           e->rt_props->debug_radiation_absorbed_this_step, */
  /*           e->rt_props->debug_star_injection_prep_iacts_with_parts_this_step, */
  /*           e->rt_props->debug_part_injection_prep_iacts_with_stars_this_step */
  /*         ); */
  /* message("Over lifetime: %12lld %12lld %12lld %12lld", */
  /*         e->rt_props->debug_radiation_emitted_tot, */
  /*         e->rt_props->debug_radiation_absorbed_tot, */
  /*         e->rt_props->debug_star_injection_prep_iacts_with_parts_tot, */
  /*         e->rt_props->debug_part_injection_prep_iacts_with_stars_tot */
  /*         ); */

  /* Have we accidentally invented or deleted some radiation somewhere? */
  if ((e->rt_props->debug_radiation_emitted_this_step !=
       e->rt_props->debug_radiation_absorbed_this_step) ||
      (e->rt_props->debug_radiation_emitted_tot !=
       e->rt_props->debug_radiation_absorbed_tot))
    error(
        "Emitted and absorbed radiation vary.\n"
        "  This step: star emission %12d; gas absorption %12d\n"
        "Since start: star emission %12lld; gas absorption %12lld",
        e->rt_props->debug_radiation_emitted_this_step,
        e->rt_props->debug_radiation_absorbed_this_step,
        e->rt_props->debug_radiation_emitted_tot,
        e->rt_props->debug_radiation_absorbed_tot);

  if ((e->rt_props->debug_part_injection_prep_iacts_with_stars_this_step !=
       e->rt_props->debug_star_injection_prep_iacts_with_parts_this_step) ||
      (e->rt_props->debug_part_injection_prep_iacts_with_stars_tot !=
       e->rt_props->debug_star_injection_prep_iacts_with_parts_tot))
    error(
        "Injection prep counts parts vs stars disagree.\n"
        "  This step: star iacts: %12d; gas iacts: %12d\n"
        "Since start: star iacts: %12lld; gas iacts: %12lld",
        e->rt_props->debug_star_injection_prep_iacts_with_parts_this_step,
        e->rt_props->debug_part_injection_prep_iacts_with_stars_this_step,
        e->rt_props->debug_star_injection_prep_iacts_with_parts_tot,
        e->rt_props->debug_part_injection_prep_iacts_with_stars_tot);

  if ((e->rt_props->debug_part_injection_prep_iacts_with_stars_this_step !=
       e->rt_props->debug_radiation_emitted_this_step) ||
      (e->rt_props->debug_part_injection_prep_iacts_with_stars_tot !=
       e->rt_props->debug_radiation_emitted_tot))
    error(
        "Injection prep iact counts vs actual iact counts disagree.\n"
        "  This step: prep iacts: %12d; inject iacts: %12d\n"
        "Since start: prep iacts: %12lld; inject iacts: %12lld",
        e->rt_props->debug_part_injection_prep_iacts_with_stars_this_step,
        e->rt_props->debug_radiation_emitted_this_step,
        e->rt_props->debug_part_injection_prep_iacts_with_stars_tot,
        e->rt_props->debug_radiation_emitted_tot);


  /* Write down energy budget for this timestep. */
  if (e->step > 1) {
    fprintf(e->rt_props->conserved_energy_filep, "\n");
  } else {
    fprintf(e->rt_props->conserved_energy_filep, "# nstars: %lld\n",
            e->total_nr_sparts);
  }
  fprintf(e->rt_props->conserved_energy_filep, "%6d %12.6e ", e->step, e->time);
  for (int g = 0; g < RT_NGROUPS; g++)
    fprintf(e->rt_props->conserved_energy_filep, "%12.6e ",
            e->rt_props->debug_total_radiation_conserved_energy[g]);

  if (e->step > 1) {
    fprintf(e->rt_props->energy_density_filep, "\n");
  } else {
    fprintf(e->rt_props->energy_density_filep, "# nstars: %lld\n",
            e->total_nr_sparts);
  }
  fprintf(e->rt_props->energy_density_filep, "%6d %12.6e ", e->step, e->time);
  for (int g = 0; g < RT_NGROUPS; g++)
    fprintf(e->rt_props->energy_density_filep, "%12.6e ",
            e->rt_props->debug_total_radiation_energy_density[g]);

  if (e->step > 1) {
    fprintf(e->rt_props->star_emitted_energy_filep, "\n");
  } else {
    fprintf(e->rt_props->star_emitted_energy_filep, "# nstars: %lld\n",
            e->total_nr_sparts);
  }
  fprintf(e->rt_props->star_emitted_energy_filep, "%6d %12.6e ", e->step,
          e->time);
  for (int g = 0; g < RT_NGROUPS; g++)
    fprintf(e->rt_props->star_emitted_energy_filep, "%12.6e ",
            e->rt_props->debug_total_star_emitted_energy[g]);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief This function is intended for debugging purposes only. It is called
 * during the self injection tasks, (regardless whether the particle actually
 * has neighbours to interact with) and intended to mark star or gas particles
 * to have been called during the step so further checks can be performed
 * further down the task system.
 *
 * @param p Hydro particle.
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_check_injection_part(struct part *restrict p,
                                  struct rt_props *props) {

  /* skip this for GEAR */
  /* if (props->do_all_parts_have_stars_checks) */
  /*   p->rt_data.debug_injection_check += 1; */
}

/**
 * @brief This function is intended for debugging purposes only. It is called
 * during the self injection tasks, (regardless whether the particle actually
 * has neighbours to interact with) and intended to mark star or gas particles
 * to have been called during the step so further checks can be performed
 * further down the task system.
 *
 * @param s Star particle.
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_check_injection_spart(struct spart *restrict s,
                                   struct rt_props *props) {

  /* skip this for GEAR */
  /* if (props->do_all_parts_have_stars_checks) */
  /*   s->rt_data.debug_injection_check += 1; */
}

/**
 * @brief Mark that a particle has been called during the gradient calls
 *
 * @param p Particle
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_count_gradient_call(struct part *restrict p) {

  p->rt_data.debug_calls_iact_gradient += 1;
}

/**
 * @brief Mark that a particle has been called during the transport calls
 *
 * @param p Particle
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_count_transport_call(struct part *restrict p) {

  p->rt_data.debug_calls_iact_transport += 1;
}
#endif /* SWIFT_RT_DEBUG_CHECKS */
#endif /* SWIFT_RT_DEBUGGING_GEAR_H */
