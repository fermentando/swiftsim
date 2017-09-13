/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_GRAVITY_IACT_H
#define SWIFT_DEFAULT_GRAVITY_IACT_H

/* Includes. */
#include "kernel_gravity.h"
#include "kernel_long_gravity.h"
#include "multipole.h"

/**
 * @brief Computes the intensity of the force at a point generated by a
 * point-mass.
 *
 * The returned quantity needs to be multiplied by the distance vector to obtain
 * the force vector.
 *
 * @param r2 Square of the distance to the point-mass.
 * @param h2 Square of the softening length.
 * @param h_inv Inverse of the softening length.
 * @param h_inv3 Cube of the inverse of the softening length.
 * @param mass Mass of the point-mass.
 * @param f_ij (return) The force intensity.
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pp_full(
    float r2, float h2, float h_inv, float h_inv3, float mass, float *f_ij) {

  /* Get the inverse distance */
  const float r_inv = 1.f / sqrtf(r2);

  /* Should we soften ? */
  if (r2 >= h2) {

    /* Get Newtonian gravity */
    *f_ij = mass * r_inv * r_inv * r_inv;

  } else {

    const float r = r2 * r_inv;
    const float ui = r * h_inv;
    float W_ij;

    kernel_grav_eval(ui, &W_ij);

    /* Get softened gravity */
    *f_ij = mass * h_inv3 * W_ij;
  }
}

/**
 * @brief Computes the intensity of the force at a point generated by a
 * point-mass truncated for long-distance periodicity.
 *
 * The returned quantity needs to be multiplied by the distance vector to obtain
 * the force vector.
 *
 * @param r2 Square of the distance to the point-mass.
 * @param h2 Square of the softening length.
 * @param h_inv Inverse of the softening length.
 * @param h_inv3 Cube of the inverse of the softening length.
 * @param mass Mass of the point-mass.
 * @param rlr_inv Inverse of the mesh smoothing scale.
 * @param f_ij (return) The force intensity.
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pp_truncated(
    float r2, float h2, float h_inv, float h_inv3, float mass, float rlr_inv,
    float *f_ij) {

  /* Get the inverse distance */
  const float r_inv = 1.f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Should we soften ? */
  if (r2 >= h2) {

    /* Get Newtonian gravity */
    *f_ij = mass * r_inv * r_inv * r_inv;

  } else {

    const float r = r2 * r_inv;
    const float ui = r * h_inv;
    float W_ij;

    kernel_grav_eval(ui, &W_ij);

    /* Get softened gravity */
    *f_ij = mass * h_inv3 * W_ij;
  }

  /* Get long-range correction */
  const float u_lr = r * rlr_inv;
  float corr_lr;
  kernel_long_grav_eval(u_lr, &corr_lr);
  *f_ij *= corr_lr;
}

/**
 * @brief Computes the force at a point generated by a multipole.
 *
 * This uses the quadrupole terms only and defaults to the monopole if
 * the code is compiled with low-order gravity only.
 *
 * @param r_x x-component of the distance vector to the multipole.
 * @param r_y y-component of the distance vector to the multipole.
 * @param r_z z-component of the distance vector to the multipole.
 * @param r2 Square of the distance vector to the multipole.
 * @param h The softening length.
 * @param h_inv Inverse of the softening length.
 * @param m The multipole.
 * @param f_x (return) The x-component of the acceleration.
 * @param f_y (return) The y-component of the acceleration.
 * @param f_z (return) The z-component of the acceleration.
 */
__attribute__((always_inline)) INLINE static void runner_iact_grav_pm(
    float r_x, float r_y, float r_z, float r2, float h, float h_inv,
    const struct multipole *m, float *f_x, float *f_y, float *f_z) {

#if SELF_GRAVITY_MULTIPOLE_ORDER < 3
  runner_iact_grav_pp_full(r2, h * h, h_inv, h_inv3, m->M_000, f_ij);
#else

  /* Get the inverse distance */
  const float r_inv = 1.f / sqrtf(r2);

  struct potential_derivatives_M2P pot;
  compute_potential_derivatives_M2P(r_x, r_y, r_z, r2, r_inv, h, h_inv, &pot);

  /* 1st order terms (monopole) */
  *f_x = m->M_000 * pot.D_100;
  *f_y = m->M_000 * pot.D_010;
  *f_z = m->M_000 * pot.D_001;

  /* 3rd order terms (quadrupole) */
  *f_x += m->M_200 * pot.D_300 + m->M_020 * pot.D_120 + m->M_002 * pot.D_102;
  *f_y += m->M_200 * pot.D_210 + m->M_020 * pot.D_030 + m->M_002 * pot.D_012;
  *f_z += m->M_200 * pot.D_201 + m->M_020 * pot.D_021 + m->M_002 * pot.D_003;
  *f_x += m->M_110 * pot.D_210 + m->M_101 * pot.D_201 + m->M_011 * pot.D_111;
  *f_y += m->M_110 * pot.D_120 + m->M_101 * pot.D_111 + m->M_011 * pot.D_021;
  *f_z += m->M_110 * pot.D_111 + m->M_101 * pot.D_102 + m->M_011 * pot.D_012;

#endif
}

#endif /* SWIFT_DEFAULT_GRAVITY_IACT_H */
