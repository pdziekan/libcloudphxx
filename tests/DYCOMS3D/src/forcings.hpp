#pragma once

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::buoyancy(const blitz::Array<real_t, 3> &th)
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;

  const real_t g = 9.81; 

//  namespace moist_air = libcloudphxx::common::moist_air;
//  const real_t eps = moist_air::R_v<real_t>() / moist_air::R_d<real_t>() - 1.;
  tmp1(ijk) = g * ((th(ijk) - th_init(ijk)) / th_init(0, 0, 0));// + eps * (rv - rv_init(ijk)));

  this->xchng_sclr(tmp1, i, j, k); 
  tmp2(i, j, k) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j, k) + tmp1(i, j, k - 1));
}

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::radiation(const blitz::Array<real_t, 3> &rv)
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;
  const auto &Tht = this->state(ix::th); 

  int nx = this->mem->grid_size[0].length(); //76
  int ny = this->mem->grid_size[1].length(); //76
  int nz = this->mem->grid_size[2].length(); //76

  // index of first cell above inversion
  blitz::thirdIndex ki;
  blitz::Array<real_t, 3> tmp(rv(ijk).shape()); // TODO: init it somewhere else
  tmp  = rv(ijk) + r_l(ijk);
  k_i(i, j) = blitz::first( tmp< setup::q_i, ki) - rv.base(blitz::thirdDim); // rv and r_l have same bases (first indices), subarrays (i.e. rv(ijk)) start with the same base as original arr (i.e. rv)!

  // calc Eqs. 5 and 6 from Ackerman et al 2009
  // TODO: z-th cell will be accounted for twice (in each integral)...
  for(int x = i.first() ; x <= i.last(); ++x) // 0..75 || 0..37 38..75 || ...
  {
    for(int y = 0 ; y < ny; ++y)
    {
      for(int z = 0 ; z < nz; ++z)
      {
        F(x, y, z) =  setup::F_0 * exp(- (nz - z) * this->dk *  blitz::sum(r_l(x, y, blitz::Range(z, nz-1))));
        F(x, y, z) += setup::F_1 * exp(- ((z) - 0) * this->dk * blitz::sum(r_l(x, y, blitz::Range(0, z))));

        if(z > k_i(x, y) )
        {
          real_t z_d = (z - k_i(x, y)) * this->dk;
          F(x, y, z) += setup::c_p * setup::rho_i * setup::D * (0.25 * pow(z_d, 4./3) + k_i(x, y) * this->dk * pow(z_d, 1./3)); 
        }
      }
    }
  }
  // Eq. 3.33 from Curry and Webster
  // calculate divergence of heat flux
  blitz::Range notopbot(1, nz-2);
  tmp1(i, j, notopbot) = (F(i, j, notopbot+1) - F(i, j, notopbot-1)) / 2./ this->dk;
  tmp1(i, j, k.last()) = (F(i, j, k.last()) - F(i, j, k.last()-1)) / this->dk;   
  tmp1(i, j, 0)        = (F(i, j, 1) - F(i, j, 0)) / this->dk;                

  // smoothing
//  this->xchng_sclr(tmp1, i, j, k);
//  F(i, j, k) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j, k) + tmp1(i, j, k - 1));
  F(i, j, k) = tmp1(i, j, k);
}

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::surf_sens()
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;
  const auto &Tht = this->state(ix::th); 

  int nx = this->mem->grid_size[0].length(); //76
  int ny = this->mem->grid_size[1].length(); //76
  int nz = this->mem->grid_size[2].length(); //76
  // exponential decrease with height
  blitz::Array<real_t, 3> hgt_fctr(nx, ny, nz);
  real_t z_0 = setup::z_rlx / si::metres;
  blitz::thirdIndex ki;
  hgt_fctr = exp(- ki * this->dk / z_0) / z_0; 
  F(i, j, k) = setup::F_sens * hgt_fctr(i, j, k);
}

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::surf_latent()
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;

  int nx = this->mem->grid_size[0].length(); //76
  int ny = this->mem->grid_size[1].length(); //76
  int nz = this->mem->grid_size[2].length(); //76
  // exponential decrease with height
  blitz::Array<real_t, 3> hgt_fctr(nx, ny, nz);
  real_t z_0 = setup::z_rlx / si::metres;
  blitz::thirdIndex ki;
  hgt_fctr = exp(- ki * this->dk / z_0) / z_0; 

  // heating[W/m^3] / rhod[kg/m^3] / latent heat of evaporation [J/kg]
  F(i, j, k) =  setup::F_lat /                           // heating 
    (libcloudphxx::common::const_cp::l_tri<real_t>() * si::kilograms / si::joules) / // latent heat of evaporation
    rhod(i, j, k) *                                                     // density
    hgt_fctr(i, j, k); 
}

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::subsidence(const int &type) // large-scale vertical wind
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;
//  tmp1(ijk) = this->state(type)(ijk);
//  this->xchng_sclr(tmp1, i, j, k);
//  tmp2(i, j, k) = - w_LS(i, j, k) * (tmp1(i, j, k + 1) - tmp1(i, j, k - 1)) / (2. * this->dk); // use built-in blitz stencil?
//  tmp1(i, j, k) = 0.25 * (tmp2(i, j, k + 1) + 2 * tmp2(i, j, k) + tmp2(i, j, k - 1));

  tmp2(ijk) = this->state(type)(ijk);
  int nz = this->mem->grid_size[2].length(); 
  blitz::Range notopbot(1, nz-2);
  tmp1(i, j, notopbot) = - w_LS(i, j, notopbot+1) * (tmp2(i, j, notopbot+1) - tmp2(i, j, notopbot-1)) / 2./ this->dk;
  tmp1(i, j, k.last()) = - w_LS(i, j, k.last()) * (tmp2(i, j, k.last()) - tmp2(i, j, k.last()-1)) / this->dk;   
  tmp1(i, j, 0)        = - w_LS(i, j, 0) * (tmp2(i, j, 1) - tmp2(i, j, 0)) / this->dk;                
}
