/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_GRADIENTS_GEAR_H
#define SWIFT_RT_GRADIENTS_GEAR_H

/* better safe than sorry */
#ifndef GIZMO_MFV_SPH
#error "Cannot compile GEAR-RT without gizmo-mfv hydro!"
#endif

#include "hydro.h" /* needed for hydro_part_geometry_well_behaved() */
#include "rt_getters.h"
#include "rt_slope_limiters_cell.h"
#include "rt_slope_limiters_face.h"

/**
 * @file src/rt/GEAR/rt_gradients.h
 * @brief Main header file for the GEAR M1 closure radiative transfer scheme
 * gradients
 */

/**
 * @brief Initialisation of the RT gradient data.
 */
__attribute__((always_inline)) INLINE static void rt_gradients_init(
    struct part *restrict p) {

  struct rt_part_data *rtd = &p->rt_data;

  for (int g = 0; g < RT_NGROUPS; g++) {
    for (int i = 0; i < 3; i++) {
      rtd->gradient[g].energy[i] = 0.f;
      rtd->gradient[g].flux[i][0] = 0.f;
      rtd->gradient[g].flux[i][1] = 0.f;
      rtd->gradient[g].flux[i][2] = 0.f;
    }
  }
}

/**
 * @brief Update the gradients for the given particle with the
 * given contributions.
 *
 * @param p Particle
 * @param g photon group index to update (0 <= g < RT_NGROUPS)
 * @param dE energy gradient contribution
 * @param dFx flux gradient contribution in x direction
 * @param dFy flux gradient contribution in y direction
 * @param dFy flux gradient contribution in z direction
 */
__attribute__((always_inline)) INLINE static void rt_gradients_update_part(
    struct part *restrict p, int g, float dE[3], float dFx[3], float dFy[3],
    float dFz[3]) {

  struct rt_part_data *rtd = &p->rt_data;

  rtd->gradient[g].energy[0] += dE[0];
  rtd->gradient[g].energy[1] += dE[1];
  rtd->gradient[g].energy[2] += dE[2];

  rtd->gradient[g].flux[0][0] += dFx[0];
  rtd->gradient[g].flux[0][1] += dFx[1];
  rtd->gradient[g].flux[0][2] += dFx[2];

  rtd->gradient[g].flux[1][0] += dFy[0];
  rtd->gradient[g].flux[1][1] += dFy[1];
  rtd->gradient[g].flux[1][2] += dFy[2];

  rtd->gradient[g].flux[2][0] += dFz[0];
  rtd->gradient[g].flux[2][1] += dFz[1];
  rtd->gradient[g].flux[2][2] += dFz[2];
}

/**
 * @brief Finalise the gradient computation after all
 * particle-particle interactions are finished.
 *
 * @param p the Particle
 **/

__attribute__((always_inline)) INLINE static void rt_finalise_gradient_part(
    struct part *restrict p) {

  /* add kernel normalization to gradients */
  const float h = p->h;
  const float h_inv = 1.0f / h;

  float norm;
  if (hydro_part_geometry_well_behaved(p)) {
    const float hinvdim = pow_dimension(h_inv);
    norm = hinvdim;
  } else {
    const float hinvdimp1 = pow_dimension_plus_one(h_inv);
    const float volume = p->geometry.volume;
    norm = hinvdimp1 * volume;
  }

  struct rt_part_data *rtd = &p->rt_data;
  for (int g = 0; g < RT_NGROUPS; g++) {
    rtd->gradient[g].energy[0] *= norm;
    rtd->gradient[g].energy[1] *= norm;
    rtd->gradient[g].energy[2] *= norm;
    for (int i = 0; i < 3; i++) {
      rtd->gradient[g].flux[i][0] *= norm;
      rtd->gradient[g].flux[i][1] *= norm;
      rtd->gradient[g].flux[i][2] *= norm;
    }
  }

/* if (p->id < 20) */
/*   message("Conserved check %.6e", rtd->conserved[0].energy); */
/*  */
/* if (p->id < 20) */
/*   message("Density check %.6e", rtd->density[0].energy); */
/*  */
/* if (p->id < 20) */
/*   message("Flux check %.6e", rtd->flux[0].energy); */
/*  */
if (p->id < 20)
  message("Gradient check %.3e %.3e %.3e | %.3e", rtd->gradient[0].energy[0], rtd->gradient[0].energy[1], rtd->gradient[0].energy[2], norm);

  rt_slope_limit_cell(p);

if (p->id < 20)
  message("Slope limit check %.6e %.6e %.6e", rtd->gradient[0].energy[0], rtd->gradient[0].energy[1], rtd->gradient[0].energy[2]);
}

/**
 * @brief symmetric gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void rt_gradients_collect(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  if (pi->rt_data.debug_injection_done != 1)
    error(
        "Trying to do symmetric iact gradient when finalise injection count is "
        "%d ID %lld",
        pi->rt_data.debug_injection_done, pi->id);

  if (pj->rt_data.debug_injection_done != 1)
    error(
        "Trying to do symmetric iact gradient when finalise injection count is "
        "%d ID %lld",
        pj->rt_data.debug_injection_done, pj->id);

  pi->rt_data.debug_calls_iact_gradient_interaction += 1;

  pj->rt_data.debug_calls_iact_gradient_interaction += 1;
#endif

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  float Bi[3][3];
  float Bj[3][3];

  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
  }

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float qi = r * hi_inv;
  kernel_deval(qi, &wi,
               &wi_dx); /* factor 1/omega for psi is swallowed in matrix */

  /* Compute kernel of pj */
  float wj, wj_dx;
  const float hj_inv = 1.0f / hj;
  const float qj = r * hj_inv;
  kernel_deval(qj, &wj,
               &wj_dx); /* factor 1/omega for psi is swallowed in matrix */

  /* Compute psi tilde */
  float psii_tilde[3];
  if (hydro_part_geometry_well_behaved(pi)) {
    psii_tilde[0] =
        wi * (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
    psii_tilde[1] =
        wi * (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
    psii_tilde[2] =
        wi * (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  } else {
    const float norm = -wi_dx * r_inv;
    psii_tilde[0] = norm * dx[0];
    psii_tilde[1] = norm * dx[1];
    psii_tilde[2] = norm * dx[2];
  }

  float psij_tilde[3];
  if (hydro_part_geometry_well_behaved(pj)) {
    psij_tilde[0] =
        wi * (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
    psij_tilde[1] =
        wi * (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
    psij_tilde[2] =
        wi * (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);
  } else {
    const float norm = -wj_dx * r_inv;
    psij_tilde[0] = norm * dx[0];
    psij_tilde[1] = norm * dx[1];
    psij_tilde[2] = norm * dx[2];
  }

  for (int g = 0; g < RT_NGROUPS; g++) {

    float Qi[4], Qj[4];
    rt_part_get_density_vector(pi, g, Qi);
    rt_part_get_density_vector(pj, g, Qj);
    const float dQ[4] = {Qi[0] - Qj[0], Qi[1] - Qj[1], Qi[2] - Qj[2],
                         Qi[3] - Qj[3]};

    /* First to the gradients of pi */
    float dE_i[3], dFx_i[3], dFy_i[3], dFz_i[3];

    /* Compute gradients for pi */
    /* there is a sign difference w.r.t. eqn. (6) because of the inverse
     * definition of dx */
    dE_i[0] = dQ[0] * psii_tilde[0];
    dE_i[1] = dQ[0] * psii_tilde[1];
    dE_i[2] = dQ[0] * psii_tilde[2];

    dFx_i[0] = dQ[1] * psii_tilde[0];
    dFx_i[1] = dQ[1] * psii_tilde[1];
    dFx_i[2] = dQ[1] * psii_tilde[2];
    dFy_i[0] = dQ[2] * psii_tilde[0];
    dFy_i[1] = dQ[2] * psii_tilde[1];
    dFy_i[2] = dQ[2] * psii_tilde[2];
    dFz_i[0] = dQ[3] * psii_tilde[0];
    dFz_i[1] = dQ[3] * psii_tilde[1];
    dFz_i[2] = dQ[3] * psii_tilde[2];

    rt_gradients_update_part(pi, g, dE_i, dFx_i, dFy_i, dFz_i);
    rt_slope_limit_cell_collect(pi, pj);

    /* Now do the gradients of pj */
    float dE_j[3], dFx_j[3], dFy_j[3], dFz_j[3];

    /* We don't need a sign change here: both the dx and the dQ
     * should switch their sign, resulting in no net change */
    dE_j[0] = dQ[0] * psij_tilde[0];
    dE_j[1] = dQ[0] * psij_tilde[1];
    dE_j[2] = dQ[0] * psij_tilde[2];

    dFx_j[0] = dQ[1] * psij_tilde[0];
    dFx_j[1] = dQ[1] * psij_tilde[1];
    dFx_j[2] = dQ[1] * psij_tilde[2];
    dFy_j[0] = dQ[2] * psij_tilde[0];
    dFy_j[1] = dQ[2] * psij_tilde[1];
    dFy_j[2] = dQ[2] * psij_tilde[2];
    dFz_j[0] = dQ[3] * psij_tilde[0];
    dFz_j[1] = dQ[3] * psij_tilde[1];
    dFz_j[2] = dQ[3] * psij_tilde[2];

    rt_gradients_update_part(pj, g, dE_j, dFx_j, dFy_j, dFz_j);
    rt_slope_limit_cell_collect(pj, pi);
  }
}

/**
 * @brief Non-symmetric gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void rt_gradients_nonsym_collect(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj) {

#ifdef SWIFT_RT_DEBUG_CHECKS
  if (pi->rt_data.debug_injection_done != 1)
    error(
        "Trying to do nonsym iact gradients when finalise injection count is "
        "%d ID %lld",
        pi->rt_data.debug_injection_done, pi->id);

  pi->rt_data.debug_calls_iact_gradient_interaction += 1;
#endif

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  float Bi[3][3];

  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
    }
  }

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float qi = r * hi_inv;
  kernel_deval(qi, &wi,
               &wi_dx); /* factor 1/omega for psi is swallowed in matrix */

  /* Compute psi tilde */
  float psii_tilde[3];
  if (hydro_part_geometry_well_behaved(pi)) {
    psii_tilde[0] =
        wi * (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
    psii_tilde[1] =
        wi * (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
    psii_tilde[2] =
        wi * (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  } else {
    const float norm = -wi_dx * r_inv;
    psii_tilde[0] = norm * dx[0];
    psii_tilde[1] = norm * dx[1];
    psii_tilde[2] = norm * dx[2];
  }

  for (int g = 0; g < RT_NGROUPS; g++) {

    float Qi[4], Qj[4];
    rt_part_get_density_vector(pi, g, Qi);
    rt_part_get_density_vector(pj, g, Qj);
    const float dQ[4] = {Qi[0] - Qj[0], Qi[1] - Qj[1], Qi[2] - Qj[2],
                         Qi[3] - Qj[3]};

    float dE_i[3], dFx_i[3], dFy_i[3], dFz_i[3];

    /* Compute gradients for pi */
    /* there is a sign difference w.r.t. eqn. (6) because of the inverse
     * definition of dx */
    dE_i[0] = dQ[0] * psii_tilde[0];
    dE_i[1] = dQ[0] * psii_tilde[1];
    dE_i[2] = dQ[0] * psii_tilde[2];

    dFx_i[0] = dQ[1] * psii_tilde[0];
    dFx_i[1] = dQ[1] * psii_tilde[1];
    dFx_i[2] = dQ[1] * psii_tilde[2];
    dFy_i[0] = dQ[2] * psii_tilde[0];
    dFy_i[1] = dQ[2] * psii_tilde[1];
    dFy_i[2] = dQ[2] * psii_tilde[2];
    dFz_i[0] = dQ[3] * psii_tilde[0];
    dFz_i[1] = dQ[3] * psii_tilde[1];
    dFz_i[2] = dQ[3] * psii_tilde[2];

    rt_gradients_update_part(pi, g, dE_i, dFx_i, dFy_i, dFz_i);

    rt_slope_limit_cell_collect(pi, pj);
  }
}

/**
 * @brief Extrapolate the given gradient over the given distance.
 *
 * @param dQ Gradient of the quantity
 * @param dx Distance vector
 * @return Change in the quantity after a displacement along the given distance
 * vector.
 */
__attribute__((always_inline)) INLINE static float rt_gradients_extrapolate(
    const float dQ[3], const float dx[3]) {

  return dQ[0] * dx[0] + dQ[1] * dx[1] + dQ[2] * dx[2];
}

/**
 * @brief Gradients reconstruction. Predict the value at point x_ij given
 * current values at particle positions and gradients at particle positions.
 */
__attribute__((always_inline)) INLINE static void rt_gradients_predict(
    struct part *restrict pi, struct part *restrict pj, float hi, float hj,
    const float *dx, float r, const float *xij_i, int group, float Qi[4],
    float Qj[4]) {

  /* Compute interface position (relative to pj, since we don't need the actual
   * position) eqn. (8)
   * Do it this way in case dx contains periodicity corrections already */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  float dE_i[3], dFx_i[3], dFy_i[3], dFz_i[3];
  float dE_j[3], dFx_j[3], dFy_j[3], dFz_j[3];
  rt_part_get_gradients(pi, group, dE_i, dFx_i, dFy_i, dFz_i);
  rt_part_get_gradients(pj, group, dE_j, dFx_j, dFy_j, dFz_j);

  float dQi[4];
  dQi[0] = rt_gradients_extrapolate(dE_i, xij_i);
  dQi[1] = rt_gradients_extrapolate(dFx_i, xij_i);
  dQi[2] = rt_gradients_extrapolate(dFy_i, xij_i);
  dQi[3] = rt_gradients_extrapolate(dFz_i, xij_i);

  float dQj[4];
  dQj[0] = rt_gradients_extrapolate(dE_j, xij_j);
  dQj[1] = rt_gradients_extrapolate(dFx_j, xij_j);
  dQj[2] = rt_gradients_extrapolate(dFy_j, xij_j);
  dQj[3] = rt_gradients_extrapolate(dFz_j, xij_j);

  /* Apply the slope limiter at this interface */
  rt_slope_limit_face(Qi, Qj, dQi, dQj, xij_i, xij_j, r);

  Qi[0] += dQi[0];
  Qi[1] += dQi[1];
  Qi[2] += dQi[2];
  Qi[3] += dQi[3];

  Qj[0] += dQj[0];
  Qj[1] += dQj[1];
  Qj[2] += dQj[2];
  Qj[3] += dQj[3];

  /* gizmo_check_physical_quantities("density", "pressure", Wi[0], Wi[1], Wi[2],
   */
  /*                                 Wi[3], Wi[4]); */
  /* gizmo_check_physical_quantities("density", "pressure", Wj[0], Wj[1], Wj[2],
   */
  /* Wj[3], Wj[4]); */
}

#endif /* SWIFT_RT_GRADIENT_GEAR_H */
