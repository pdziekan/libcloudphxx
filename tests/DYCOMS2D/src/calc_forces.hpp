#pragma once
#include "forcings.hpp"

template <class ct_params_t>
void kin_cloud_2d_lgrngn<ct_params_t>::rv_src()
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  int nz = this->mem->grid_size[1].length(); //76
  // surface flux
  surf_latent();
  // divergence of rv flux
  blitz::Range notop(0, nz-2);
  alpha(i, notop) = (F(i, notop) - F(i, notop+1)) / this->dj;
  alpha(i, j.last()) = F(i, j.last()) / this->dj;
  // change of rv[1/s] = latent heating[W/m^3] / lat_heat_of_evap[J/kg] / density[kg/m^3]
  alpha(ijk)/=(libcloudphxx::common::const_cp::l_tri<real_t>() * si::kilograms / si::joules) * rhod(ijk);

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
  // divergence of th flux, F(j) is upward flux through bottom of the j-th cell
  int nz = this->mem->grid_size[1].length(); //76
  blitz::Range notop(0, nz-2);
  alpha(i, notop) = (F(i, notop) - F(i, notop+1)) / this->dj;
  alpha(i, j.last()) = F(i, j.last()) / this->dj;
  // radiation
  radiation(rv);
  // divergence of th flux, F(j) is upward flux in the middle of the j-th cell
  blitz::Range notopbot(1, nz-2);
  alpha(i, notopbot) += ( -F(i, notopbot+1) + F(i, notopbot-1)) / 2./ this->dj;
  alpha(i, j.last()) += ( -F(i, j.last()) + F(i, j.last()-1)) / this->dj;
  alpha(i, 0)        += ( -F(i, 1) + F(i, 0)) / this->dj;

  // change of theta[K/s] = heating[W/m^3] * theta[K] / T[K] / c_p[J/K/kg] / rhod[kg/m^3]
  for(int x = i.first() ; x <= i.last(); ++x)
  {
      for(int z = j.first() ; z <= j.last(); ++z)
      {
        alpha(x, z) = alpha(x, z) * this->state(ix::th)(x, z) / rhod(x, z) /
                     (libcloudphxx::common::moist_air::c_p<real_t>(this->state(ix::rv)(x, z)) * si::kilograms * si::kelvins / si::joules) /
                     (libcloudphxx::common::theta_dry::T<real_t>(this->state(ix::th)(x, z) * si::kelvins, rhod(x, z) * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins);
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
void kin_cloud_2d_lgrngn<ct_params_t>::w_src(const blitz::Array<real_t, 2> &th)
{
  const auto &ijk = this->ijk;
  // buoyancy
  buoyancy(th);
  alpha(ijk) = F(ijk);
  // large-scale vertical wind
  subsidence(ix::w);
  alpha(ijk) += F(ijk);
}


