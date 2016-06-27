#pragma once
#include "forcings.hpp"

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::rv_src()
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;

  if(this->timestep >= this->spinup)
  {
    // surface flux
    surf_latent();
    // sum of rv flux
    int nz = this->mem->grid_size[2].length();
    blitz::Range notop(0, nz-2);
    alpha(i, j, notop) = (F(i, j, notop) - F(i, j, notop+1)) / this->dk;
    alpha(i, j, k.last()) = F(i, j, k.last()) / this->dk;

    // top and bottom cells are two times lower
    alpha(i, j, 0) *= 2; 
    alpha(i, j, this->k.last()) *= 2; 
    // change of rv[1/s] = latent heating[W/m^3] / lat_heat_of_evap[J/kg] / density[kg/m^3]
    alpha(ijk)/=(libcloudphxx::common::const_cp::l_tri<real_t>() * si::kilograms / si::joules) * this->rhod(ijk);
  }
  else
    alpha(ijk) = 0;

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
  // beta as tmp storage
  beta(ijk) = F(ijk);
  // radiation
  radiation(rv);
  // add fluxes from radiation and surface
  F(ijk) += beta(ijk);

  // sum of th flux, F(j) is upward flux through the bottom of the j-th cell

  int nz = this->mem->grid_size[2].length(); //76
  // sum of th flux
  blitz::Range notop(0, nz-2);
  alpha(i, j, notop) = (F(i, j, notop) - F(i, j, notop+1)) / this->dk;
  alpha(i, j, k.last()) = F(i, j, k.last()) / this->dk;

  // top and bottom cells are two times lower
  alpha(i, j, 0) *= 2; 
  alpha(i, j, this->j.last()) *= 2; 

  // change of theta[K/s] = heating[W/m^3] * theta[K] / T[K] / c_p[J/K/kg] / rhod[kg/m^3]
  for(int x = i.first() ; x <= i.last(); ++x)
  {
    for(int y = j.first() ; y <= j.last(); ++y)
    {
      for(int z = k.first() ; z <= k.last(); ++z)
      {
        alpha(x, y, z) = alpha(x, y, z) * this->state(ix::th)(x, y, z) / this->rhod(x, y, z) /
                     (libcloudphxx::common::moist_air::c_p<real_t>(rv(x, y, z)) * si::kilograms * si::kelvins / si::joules) /
                     (libcloudphxx::common::theta_dry::T<real_t>(this->state(ix::th)(x, y, z) * si::kelvins, this->rhod(x, y, z) * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins);
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
void kin_cloud_3d_lgrngn<ct_params_t>::w_src(const blitz::Array<real_t, 3> &th, const blitz::Array<real_t, 3> &rv)
{
  const auto &ijk = this->ijk;
  // buoyancy
  buoyancy(th, rv);
  alpha(ijk) = F(ijk);
  // large-scale vertical wind
  subsidence(ix::w);
  alpha(ijk) += F(ijk);
}


