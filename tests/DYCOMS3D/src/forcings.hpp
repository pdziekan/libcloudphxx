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
  F(i, j, k) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j, k) + tmp1(i, j, k - 1));
}

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::radiation(const blitz::Array<real_t, 3> &rv)
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;

  int ny = this->mem->grid_size[1].length(); //76
  int nz = this->mem->grid_size[2].length(); //76

  // index of first cell above inversion
  blitz::thirdIndex ki;
  tmp1(ijk)  = rv(ijk) + r_l(ijk);
  k_i(i, j) = blitz::first( tmp1 < setup::q_i, ki); 
  F(ijk) = 0;

  // calc Eqs. 5 and 6 from Ackerman et al 2009
  // TODO: z-th cell will be accounted for twice (in each integral)...
  for(int x = i.first() ; x <= i.last(); ++x) // 0..75 || 0..37 38..75 || ...
  {
    for(int y = 0 ; y < ny; ++y)
    {
      for(int z = 0 ; z < nz; ++z)
      {
        real_t sum = r_l(x, y, z) / 2.;
        if(z < nz-1) sum += blitz::sum(r_l(x, y, blitz::Range(z+1, nz-1)));
        F(x, y, z) += setup::F_0 * exp(- (nz - z) * this->dk * sum); 

        sum = r_l(x, y, z) / 2.;
        if(z > 0) sum += blitz::sum(r_l(x, y, blitz::Range(0, z-1)));
        F(x, y, z) += setup::F_1 * exp(- ((z) - 0) * this->dk * sum);

        if(z >= k_i(x, y) )
        {
          real_t z_i = (-0.5 + k_i(x,y)) * this->dk;
          real_t z_d = z * this->dk - z_i;
          F(x, y, z) += setup::c_p * setup::rho_i * setup::D * (0.25 * pow(z_d, 4./3) + z_i * pow(z_d, 1./3)); 
        }
      }
    }
  }
}

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::surf_sens()
{
  const auto &ijk = this->ijk;
  F(ijk) = setup::F_sens * hgt_fctr(ijk);
}

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::surf_latent()
{
  const auto &ijk = this->ijk;
  F(ijk) =  setup::F_lat * hgt_fctr(ijk); 
}

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::subsidence(const int &type) // large-scale vertical wind
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;
  tmp1(ijk) = this->state(type)(ijk);
  this->xchng_sclr(tmp1, i, j, k);
  F(i, j, k) = - w_LS(i, j, k) * (tmp1(i, j, k + 1) - tmp1(i, j, k - 1)) / (2. * this->dk); // use built-in blitz stencil?
}