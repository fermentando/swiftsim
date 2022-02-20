/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c)  2022 Tsang Keung Chan (chantsangkeung@gmail.com)
 *                2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_UNPHYSICAL_SPHM1RT_H
#define SWIFT_RT_UNPHYSICAL_SPHM1RT_H

/**
 * @file src/rt/SPHM1RT/rt_unphysical.h
 * @brief Routines for checking for and correcting unphysical scenarios
 */

#include <string.h>

#include "rt_struct.h"
#include "rt_properties.h"

/**
 * @brief check for and correct if needed unphysical
 * values for a radiation state.
 *
 * @param energy_density pointer to the radiation energy density
 * @param flux pointer to radiation flux (3 dimensional)
 * @param e_old energy density before change to check. Set = 0 if not available
 * @param cred reduced speed of light in code unit
 */
__attribute__((always_inline)) INLINE static void rt_check_unphysical_state(
    float* energy_density, float* flux, const float e_old, const float cred) {

  /* Check for negative energies */
  /* Note to self for printouts: Maximal allowable F = E * c.
   * In some cases, e.g. while cooling, we don't modify the fluxes,
   * so you can get an estimate of what the photon energy used to be
   * by dividing the printed out fluxes by the speed of light in
   * code units */
  if (isinf(*energy_density) || isnan(*energy_density))
    error("Got inf/nan radiation energy case | %.6e | %.6e %.6e %.6e",
          *energy_density, flux[0], flux[1], flux[2]);

  if (*energy_density <= 0.f) {
    *energy_density = 0.f;
    flux[0] = 0.f;
    flux[1] = 0.f;
    flux[2] = 0.f;
    return;
  }

  /* Check for too high fluxes */
  const float flux2 = flux[0] * flux[0] + flux[1] * flux[1] + flux[2] * flux[2];

  if (isinf(flux2) || isnan(flux2))
    error("Got inf/nan in flux2 | %.6e| %.6e %.6e %.6e", flux2, flux[0],
          flux[1], flux[2]);

  const float flux_norm = (flux2 == 0.f) ? 0.f : sqrtf(flux2);
  const float flux_norm_inv = (flux_norm == 0.f) ? 0.f : 1.f / flux_norm;
  const float flux_max = cred * *energy_density;
  float flux_diff = flux_norm - flux_max;

  if (flux_norm != 0.f) {
    if (flux_diff > 0.f) {
      const float correct = flux_max * flux_norm_inv;
      flux[0] *= correct;
      flux[1] *= correct;
      flux[2] *= correct;
    }
  }
}


/**
 * @brief check whether gas species abundances have physical
 * values and correct small errors if necessary.
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void
rt_check_unphysical_abundances(struct part* restrict p) {



  double test_abundance;
  char name[10];
  for (int spec = 0; spec < rt_species_count; spec++) {
    test_abundance = p->rt_data.tchem.abundances[spec];
    if (test_abundance < 0.f) {
      if (test_abundance < -1e4) {
        sprintf(name, "%s",
            rt_cooling_get_species_name((enum rt_cooling_species)spec));
        message("WARNING: Got negative abundance in %s", name);
      }
      p->rt_data.tchem.abundances[spec] = 0.f;
    }    
  }

  float abundance_tot = 0.f;
  for (int j = 0; j < rt_species_count; j++) {
    abundance_tot += p->rt_data.tchem.abundances[j];
  }

  /* Make sure we sum up to 1. TODO: Assuming we have no metals. */
  if (fabsf(abundance_tot - 1.f) > 1e-3)
    error("Got total mass fraction of gas = %.6g", abundance_tot);
}

#endif /* SWIFT_RT_UNPHYSICAL_SPHM1RT_H */
