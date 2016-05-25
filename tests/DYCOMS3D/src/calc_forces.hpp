#pragma once
#include "forcings.hpp"

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::rv_src()
{
  const auto &ijk = this->ijk;
  // surface flux
  surf_latent();
  alpha(ijk) = F(ijk);
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
void kin_cloud_3d_lgrngn<ct_params_t>::th_src(const blitz::Array<real_t, 3> &rv)
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;
  // -- heating --
  // surface flux
  surf_sens();
  alpha(ijk) = F(ijk);
  // radiation
  radiation(rv);
  alpha(ijk) += F(ijk);
  // change of theta[K/s] = heating[W/m^3] * theta[K] / T[K] / c_p[J/K/kg] / rhod[kg/m^3]
  for(int x = i.first() ; x <= i.last(); ++x)
  {
    for(int y = j.first() ; y <= j.last(); ++y)
    {
      for(int z = k.first() ; z <= k.last(); ++z)
      {
        alpha(x, y, z) = alpha(x, y, z) * this->state(ix::th)(x, y, z) / rhod(x, y, z) /
                     (libcloudphxx::common::moist_air::c_p<real_t>(rv(x, y, z)) * si::kilograms * si::kelvins / si::joules) /
                     (libcloudphxx::common::theta_dry::T<real_t>(this->state(ix::th)(x, y, z) * si::kelvins, rhod(x, y, z) * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins);
      }
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
void kin_cloud_3d_lgrngn<ct_params_t>::w_src(const blitz::Array<real_t, 3> &th)
{
  const auto &ijk = this->ijk;
  // buoyancy
  buoyancy(th);
  alpha(ijk) = F(ijk);
  // large-scale vertical wind
  subsidence(ix::w);
  alpha(ijk) += F(ijk);
}


