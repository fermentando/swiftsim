/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_KIARA_CHEMISTRY_IACT_H
#define SWIFT_KIARA_CHEMISTRY_IACT_H

/* Local includes */
#include "random.h"
#include "tools.h"
#include "physical_constants.h"
#include "units.h"
#include "cooling.h"
#include "grackle.h"
#include "hydro.h"
#include "feedback.h"

/**
 * @file KIARA/chemistry_iact.h
 * @brief Smooth metal interaction functions following the KIARA model.
 */


/**
 * @brief Firehose feedback interaction between wind - gas particles (symmetric) during decoupled time.
 * Used for updating properties of gas particles neighbouring a star particle                                                                                                 
 *                                                                                                                                                                             
 * @param r2 Comoving square distance between the two particles.                                                                                                               
 * @param dx Comoving vector separating both particles (pi - pj).                                                                                                              
 * @param hi Comoving smoothing-length of particle i.                                                                                                                          
 * @param hj Comoving smoothing-length of particle j.                                                                                                                          
 * @param pi Wind particle (not updated).                                                                                                                          
 * @param pj Gas particle. 
 * @param xpi Extra particle data (wind)                                                                                                                                           
 * @param xpj Extra particle data (gas)                                                                                                                                          
 * @param cosmo The cosmological model.                                                                                                                                        
 * @param fb_props Properties of the feedback scheme.                                                                                                                          
 * @param ti_current Current integer time used value for seeding random number generator   
 * 
 * 
 * @param phys_const Physical constants
 * @param us Unit system
 * @param cooling Cooling structure                                                                                     
 * */
__attribute__((always_inline)) INLINE static void
firehose_wind(
    const float r2, const float hi, const float hj,
    struct part *pi, struct xpart *xpi, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props,
    const integertime_t ti_current,

    
    const struct phys_const* phys_const,
    const struct unit_system* us, const struct cooling_function_data* cooling, FILE *fp) {


  /* Ignore COUPLED particles */
  if (pi->feedback_data.decoupling_delay_time <= 0.f) return;
  
  /* No wind-wind interaction */
  if (pj->feedback_data.decoupling_delay_time >= 0.f) return;

  message("FIREHOSE MODEL ACTIVATED");
  //error("Reached first wind particle");
  
  /* Mixed wind particle*/
  struct part *pk = NULL;
  struct xpart *xpk = NULL;

  //bzero(pk, sizeof(struct part));
  memmove(pk, pi, sizeof(struct part));


  /* Gas particle density */
  const float rho_j = hydro_get_comoving_density(pj);

  /* Get r. */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);
  warning("FIREHOSE test kernel: %f\n", wi);

  /* Constants */ 
  const float gamma = 5/3.;
  float dt =
      hydro_compute_timestep(pi, xpi, hydro_props, cosmo); // or get_timestep(p->time_bin, e->time_base) from feedback.c

  /* Kinetic and thermal properties of interacting gas */
  const float T_threshold = 1e4 / us->UnitTemperature_in_cgs;
  if (cooling_convert_u_to_temp(pi->u, xpi->cooling_data.e_frac, cooling, pi) < T_threshold){
    pi->u = cooling_convert_temp_to_u(T_threshold, xpi->cooling_data.e_frac, cooling, pi);
  }
  float Thot = cooling_convert_u_to_temp(pj->u, xpj->cooling_data.e_frac, cooling, pj);
  float Tstream = cooling_convert_u_to_temp(pi->u, xpi->cooling_data.e_frac, cooling, pi);

  const float X_Hi = chemistry_get_metal_mass_fraction_for_cooling(pi)[chemistry_element_H];
  const float X_Hj = chemistry_get_metal_mass_fraction_for_cooling(pj)[chemistry_element_H];
  const float yhelium_i = (1. - X_Hi) / (4. * X_Hi);
  const float yhelium_j = (1. - X_Hj) / (4. * X_Hj);
  const float mu_i = (1. + yhelium_i) / (1. + xpi->cooling_data.e_frac + 4. * yhelium_i);
  const float mu_j = (1. + yhelium_j) / (1. + xpj->cooling_data.e_frac + 4. * yhelium_j);
  float chi = Thot/mu_j/Tstream*mu_i; //assuming collisional equilibrium
  

  /* We compute the velocity of the phases*/
  float stream_prior_v2 = 0;
  float surroundings_prior_v2 = 0;

  for (int i=0; i<3; i++){
    stream_prior_v2 += pi->v_full[i]*pi->v_full[i];
    surroundings_prior_v2 += pj->v_full[i]*pj->v_full[i];
  }

  float Mach = sqrt(stream_prior_v2-surroundings_prior_v2)/(sqrt(pj->u/rho_j/gamma/(gamma - 1)));  


  /* Estimate mixing layer properties from ghost particles */
  pk->rho = sqrt(pj->rho*pi->rho);    
  pk->u = cooling_convert_temp_to_u(sqrt(Tstream*Thot), xpk->cooling_data.e_frac, cooling, pk);

  const float Lambda_mix = sqrt(pk->u)
    /cooling_time(phys_const, us, hydro_props, cosmo, cooling, pk, xpk);


  /* If wind has just been ejected, start destruction time variable */

  if (pi->chemistry_data.destruction_time <= 0.f){
    pi->chemistry_data.destruction_time = 0.f;
  }

  /* Define cooling and destruction timescales for streams*/

  float tcoolmix = phys_const->const_boltzmann_k* sqrt(Tstream*Thot) / ((gamma -1 )* Lambda_mix);
  float tsc = 2 * pi->chemistry_data.radius_stream / sqrt(pi->u/chi/rho_j*gamma/(gamma - 1));
  float tshear = pi->chemistry_data.radius_stream/ sqrt(stream_prior_v2 - surroundings_prior_v2);

  float virtual_mass = chi*rho_j*pow(pi->chemistry_data.radius_stream,2)*M_PI;
  float mdot;
  float delta_m;

  /* If the cloud is destroyed, updated destruction time and mdot*/
  if (tcoolmix/tshear > 1.f){

    /* If cloud just began to be destroyed, update initial mass */
    if (pi->chemistry_data.destruction_time == 0.f || Mach < 1){
      pi->chemistry_data.initial_mass = virtual_mass;
    }

    pi->chemistry_data.destruction_time += dt;
    mdot = -1/tshear*pi->chemistry_data.initial_mass*exp(-pi->chemistry_data.destruction_time/tshear);
    delta_m = fabs(mdot);
  }

  /* If cloud survives, cancel destruction time and update mdot*/
  if (tcoolmix/tshear < 1.f) {

    /* If cloud just began to grow, update initial mass */
    if (pi->chemistry_data.destruction_time > 0.f || Mach < 1){
      pi->chemistry_data.initial_mass = virtual_mass;
    }

    pi->chemistry_data.destruction_time = 0.f;
    mdot = 4/chi * pi->chemistry_data.initial_mass/tsc * pow(tcoolmix/tsc, -0.25);
    delta_m = fabs(mdot);
  
  } 

  /* Change in properties of wind and surroundings as it travels: 

   * 1) Update chemistry */
  for (int elem = 0; elem < chemistry_element_count; ++elem) { 
    pi->chemistry_data.metal_mass_fraction[elem]  = wi* 1/pi->mass * ((pi->mass - delta_m)* pi->chemistry_data.metal_mass_fraction[elem] + delta_m * pj->chemistry_data.metal_mass_fraction[elem]);
    pj->chemistry_data.metal_mass_fraction[elem]  = wi * 1/pj->mass * ((pj->mass - delta_m)* pj->chemistry_data.metal_mass_fraction[elem] + delta_m * pi->chemistry_data.metal_mass_fraction[elem]);

  }

   /* 2) Update particles' internal energy */
   pi->u = wi * 1/pi->mass * ((pi->mass - delta_m)*pi->u + delta_m * pj->u);
   pj->u = wi * 1/pj->mass * ((pj->mass - delta_m)*pj->u + delta_m * pi->u);


   /* 3) Conserve momentum */
  float stream_post_v2 = 0;
  float surroundings_post_v2 = 0;
  for (int i=0; i<3; i++){

    pi->v_full[i] = wi * 1/pi->mass * ((pi->mass - delta_m)*pi->v_full[i] + delta_m * pj->v_full[i]);
    pj->v_full[i] = wi * 1/pj->mass * ((pj->mass - delta_m)*pj->v_full[i] + delta_m * pi->v_full[i]);

    stream_post_v2 += pi->v_full[i]*pi->v_full[i];
    surroundings_post_v2 += pj->v_full[i]*pj->v_full[i];
  }

   /* 4) Deposit excess energy onto stream */
  float delE = 0.5*(pi->mass * (stream_post_v2 - stream_prior_v2) + pj->mass*(surroundings_post_v2 - surroundings_prior_v2));
  pi->u += wi*delE/pi->mass;



  /* Update virtual mass of stream */
  virtual_mass += wi*mdot*dt;
  pi->chemistry_data.radius_stream = sqrt(virtual_mass/M_PI/chi/rho_j);

    /* Dust destruction */
  float rho_dust = pi->cooling_data.dust_mass * kernel_root *hi_inv;
  float tsp = cooling->dust_grainsize * 3.2e-18 / (pow(us->UnitLength_in_cgs, 4) / us->UnitTime_in_cgs) *
    rho_dust/phys_const->const_proton_mass / (pow(2e6/us->UnitTemperature_in_cgs/Tstream, 2.5) +1);
  
  float delta_rho = -3*rho_dust/tsp*dt;
  pi->cooling_data.dust_mass += delta_rho / (kernel_root *hi_inv);

  /* Recouple if Mach < 1 */
  if (Mach < 1 )pi->feedback_data.decoupling_delay_time = 0.f;

  if (pi->id % 1000 != 0) return;
  /* Print wind properties*/
  const float length_convert = cosmo->a * fb_props->length_to_kpc;
  const float velocity_convert = cosmo->a_inv / fb_props->kms_to_internal;
  const float rho_convert = cosmo->a3_inv * fb_props->rho_to_n_cgs;
  const float u_convert =
      cosmo->a_factor_internal_energy / fb_props->temp_to_u_factor;
  fprintf(fp, "%.3f %lld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d\n",
        cosmo->z,
        pi->id,
        pi->gpart->fof_data.group_mass * fb_props->mass_to_solar_mass, 
        pi->h * cosmo->a * fb_props->length_to_kpc,
        pi->x[0] * length_convert,
        pi->x[1] * length_convert,
        pi->x[2] * length_convert,
        pi->v_full[0] * velocity_convert,
        pi->v_full[1] * velocity_convert,
        pi->v_full[2] * velocity_convert,
        pi->u * u_convert,
        pi->rho * rho_convert,
        pi->feedback_data.radius_stream * length_convert,
        pi->chemistry_data.metal_mass_fraction_total,
        pi->viscosity.v_sig * velocity_convert,
        pi->feedback_data.decoupling_delay_time * fb_props->time_to_Myr,
        pi->feedback_data.number_of_times_decoupled);

  return;
  
}


/**
 * @brief do chemistry computation after the runner_iact_density (symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_chemistry(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  float wi, wi_dx;
  float wj, wj_dx;

  /* Get the masses. */
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute the kernel function for pj */
  const float uj = r / hj;
  kernel_deval(uj, &wj, &wj_dx);

  const float wi_dr = wi_dx * r_inv;
  const float mj_wi_dr = mj * wi_dr;

  const float wj_dr = wj_dx * r_inv;
  const float mi_wj_dr = mi * wj_dr;

    /* Compute the shear tensor */
  for (int i = 0; i < 3; i++) {
    const float dxi_mj_wi_dr = dx[i] * mj_wi_dr;
    const float dxi_mi_wj_dr = dx[i] * mi_wj_dr;

    chi->shear_tensor[i][0] += (pj->v[0] - pi->v[0]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][1] += (pj->v[1] - pi->v[1]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][2] += (pj->v[2] - pi->v[2]) * dxi_mj_wi_dr;

    chj->shear_tensor[i][0] -= (pj->v[0] - pi->v[0]) * dxi_mi_wj_dr;
    chj->shear_tensor[i][1] -= (pj->v[1] - pi->v[1]) * dxi_mi_wj_dr;
    chj->shear_tensor[i][2] -= (pj->v[2] - pi->v[2]) * dxi_mi_wj_dr;
  }
}

/**
 * @brief do chemistry computation after the runner_iact_density (non symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_chemistry(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

  struct chemistry_part_data *chi = &pi->chemistry_data;

  float wi, wi_dx;

  /* Get the masses. */
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);

  const float wi_dr = wi_dx * r_inv;
  const float mj_wi_dr = mj * wi_dr;

  /* Compute the shear tensor */
  for (int i = 0; i < 3; i++) {
    const float dxi_mj_wi_dr = dx[i] * mj_wi_dr;
    chi->shear_tensor[i][0] += (pj->v[0] - pi->v[0]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][1] += (pj->v[1] - pi->v[1]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][2] += (pj->v[2] - pi->v[2]) * dxi_mj_wi_dr;
  }
}

/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param ti_current The current time (in integer)
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *pi, struct xpart *xpi, struct part *pj, struct xpart *xpj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct chemistry_global_data* cd,

    const struct feedback_props *fb_props, const struct hydro_props *hydro_props, const struct phys_const* phys_const,
    const struct unit_system* us, const struct cooling_function_data* cooling, FILE *fp) {
     
    struct chemistry_part_data *chi = &pi->chemistry_data; 
    struct chemistry_part_data *chj = &pj->chemistry_data;

  /* No need to diffuse if both particles are not diffusing. */
  if (chj->diffusion_coefficient > 0 && chi->diffusion_coefficient > 0) {

    /* Get mass */
    const float mj = hydro_get_mass(pj);
    const float mi = hydro_get_mass(pi);
    const float rhoj = hydro_get_physical_density(pj, cosmo);
    const float rhoi = hydro_get_physical_density(pi, cosmo);

    float wi, wj, dwi_dx, dwj_dx;

    /* Get r */
    const float r = sqrtf(r2);

    /* part j */
    /* Get the kernel for hj */
    const float hj_inv = 1.0f / hj;

    /* Compute the kernel function for pj */
    const float xj = r * hj_inv;
    kernel_deval(xj, &wj, &dwj_dx);

    /* part i */
    /* Get the kernel for hi */
    const float hi_inv = 1.0f / hi;

    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);

    /* Get 1/r */
    const float r_inv = 1.f / sqrtf(r2);

    const float wi_dr = dwi_dx * r_inv;
    const float wj_dr = dwj_dx * r_inv;

    const float mj_dw_r = mj * wi_dr;
    const float mi_dw_r = mi * wj_dr;

    const float rhoij_inv = 1.f / (rhoi * rhoj);

    /**
     * Compute the diffusion following Eq. 2.14
     * from Monaghan, Huppert, & Worster (2006).
     */
    float coef = 4.f * chi->diffusion_coefficient * chj->diffusion_coefficient;
    coef /= chi->diffusion_coefficient + chj->diffusion_coefficient;

    const float coef_i = coef * mj_dw_r * rhoij_inv;
    const float coef_j = coef * mi_dw_r * rhoij_inv;

    /* Compute the time derivative of metals due to diffusion */
    const float dZ_ij_tot = chi->metal_mass_fraction_total - chj->metal_mass_fraction_total;
    chi->dZ_dt_total += coef_i * dZ_ij_tot;
    chj->dZ_dt_total -= coef_j * dZ_ij_tot;

    for (int elem = 0; elem < chemistry_element_count; elem++) {
      const float dZ_ij = chi->metal_mass_fraction[elem] - chj->metal_mass_fraction[elem];
      chi->dZ_dt[elem] += coef_i * dZ_ij;
      chj->dZ_dt[elem] -= coef_j * dZ_ij;
    }
  }
  if (cd->firehose_feedback!=0.f){
      firehose_wind(r2, hi, hj, pi, xpi, pj, xpj, cosmo, hydro_props, fb_props, t_current,
    phys_const, us, cooling, fp);
  }


}

/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (nonsymmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param ti_current The current time (in integer)
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology) {

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  if (chj->diffusion_coefficient > 0 && chi->diffusion_coefficient > 0) {

    /* Get mass */
    const float mj = hydro_get_mass(pj);
    const float rhoj = hydro_get_physical_density(pj, cosmo);
    const float rhoi = hydro_get_physical_density(pi, cosmo);

    float wi, dwi_dx;

    /* Get r */
    const float r = sqrtf(r2);

    /* part i */
    /* Get the kernel for hi */
    const float hi_inv = 1.0f / hi;

    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);

    /* Get 1/r */
    const float r_inv = 1.f / sqrtf(r2);
    const float wi_dr = dwi_dx * r_inv;

    const float mj_dw_r = mj * wi_dr;

    const float rhoij_inv = 1.f / (rhoi * rhoj);

    /**
     * Compute the diffusion following Eq. 2.14
     * from Monaghan, Huppert, & Worster (2006).
     */
    float coef = 4.f * chi->diffusion_coefficient * chj->diffusion_coefficient;
    coef /= chi->diffusion_coefficient + chj->diffusion_coefficient;

    const float coef_i = coef * mj_dw_r * rhoij_inv;

    /* Compute the time derivative */
    const float dZ_ij_tot = chi->metal_mass_fraction_total - chj->metal_mass_fraction_total;
    chi->dZ_dt_total += coef_i * dZ_ij_tot;

    for (int elem = 0; elem < chemistry_element_count; elem++) {
      const float dZ_ij = chi->metal_mass_fraction[elem] - chj->metal_mass_fraction[elem];
      chi->dZ_dt[elem] += coef_i * dZ_ij;
    }
  }
}


#endif /* SWIFT_KIARA_CHEMISTRY_IACT_H */
