#pragma once
#include "forcings.hpp"

template <class ct_params_t>
void kin_cloud_2d_lgrngn<ct_params_t>::rv_src()
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  // surface flux
  surf_latent();
  // sum of rv flux
  this->xchng_sclr(F, i, j);
  alpha(i, j) = ( F(i, j) - F(i, j+1)) / this->dj;
  // top and bottom cells are two times lower
  alpha(i, 0) *= 2; 
  alpha(i, this->j.last()) *= 2; 
  // change of rv[1/s] = latent heating[W/m^3] / lat_heat_of_evap[J/kg] / density[kg/m^3]
  alpha(ijk)/=(libcloudphxx::common::const_cp::l_tri<real_t>() * si::kilograms / si::joules) * this->rhod(ijk);
  // large-scale vertical wind
  subsidence(ix::rv);
  alpha(ijk) += F(ijk);
  // absorber
  alpha(ijk) += (*this->mem->vab_coeff)(ijk) * this->rv_eq(ijk);
  // TODO: add nudging to alpha
  beta(ijk) = - (*this->mem->vab_coeff)(ijk);
  // TODO: add nudging to beta
}


template <class ct_params_t>
void kin_cloud_2d_lgrngn<ct_params_t>::th_src(const blitz::Array<real_t, 2> &rv)
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  // -- heating --
  // surface flux
  surf_sens();
  // beta as tmp storage
  beta(ijk) = F(ijk);
  // radiation
  radiation(rv);
  // add fluxes from radiation and surface
  F(ijk) += beta(ijk);
  // sum of th flux, F(j) is upward flux through the bottom of the j-th cell
  this->xchng_sclr(F, i, j);
  alpha(i, j) = ( F(i, j) - F(i, j+1)) /  this->dj;
  // top and bottom cells are two times lower
  alpha(i, 0) *= 2; 
  alpha(i, this->j.last()) *= 2; 

  // change of theta[K/s] = heating[W/m^3] * theta[K] / T[K] / c_p[J/K/kg] / this->rhod[kg/m^3]
  for(int x = i.first() ; x <= i.last(); ++x)
  {
      for(int z = j.first() ; z <= j.last(); ++z)
      {
        alpha(x, z) = alpha(x, z) * this->state(ix::th)(x, z) / this->rhod(x, z) /
                     (libcloudphxx::common::moist_air::c_p<real_t>(rv(x, z)) * si::kilograms * si::kelvins / si::joules) /
                     (libcloudphxx::common::theta_dry::T<real_t>(this->state(ix::th)(x, z) * si::kelvins, this->rhod(x, z) * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins);
      }
  }
  // large-scale vertical wind
  subsidence(ix::th);
  alpha(ijk) += F(ijk);
  // absorber
  alpha(ijk) += (*this->mem->vab_coeff)(ijk) * this->th_eq(ijk);
  // TODO: add nudging to alpha
  beta(ijk) = - (*this->mem->vab_coeff)(ijk);
  // TODO: add nudging to beta
}

template <class ct_params_t>
void kin_cloud_2d_lgrngn<ct_params_t>::w_src(const blitz::Array<real_t, 2> &th, const blitz::Array<real_t, 2> &rv)
{
  const auto &ijk = this->ijk;
  // buoyancy
  buoyancy(th, rv);
  alpha(ijk) = F(ijk);
  // large-scale vertical wind
  subsidence(ix::w);
  alpha(ijk) += F(ijk);
}


