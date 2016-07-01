#pragma once

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::buoyancy(const blitz::Array<real_t, 3> &th, const blitz::Array<real_t, 3> &rv)
{
  const auto &ijk = this->ijk;
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;

  const real_t g = 9.81; 

  namespace moist_air = libcloudphxx::common::moist_air;
  const real_t eps = moist_air::R_v<real_t>() / moist_air::R_d<real_t>() - 1.;
  tmp1(ijk) = g * ((th(ijk) - this->th_eq(ijk)) / this->th_ref(ijk)) + eps * (rv(ijk) - this->rv_eq(ijk)
) - r_l(ijk);

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
        setup::real_t sum = blitz::sum(r_l(x, y, blitz::Range(z, nz-1)));
        if(z==0)
          F(x, y, z) += setup::F_0 * exp(- (nz - z - 1) * this->dk * sum); 
        else
          F(x, y, z) += setup::F_0 * exp(- (nz - z - 0.5) * this->dk * sum); 
   
        if(z > 0)
        {
          sum = blitz::sum(r_l(x, y, blitz::Range(0, z-1)));
          F(x, y, z) += setup::F_1 * exp(- (z - 0.5) * this->dk * sum);
        }
  
        if(z > k_i(x, y) )
        {
          real_t z_i = (k_i(x, y) - .5) * this->dk; // bottom of the first cell above inversion, z=0 at k=0.5
          real_t z_d = (z - 0.5) * this->dk - z_i;
          F(x, y, z) += setup::c_p * setup::rho_i * setup::D * (0.25 * pow(z_d, 4./3) + z_i * pow(z_d, 1./3)); 
        }
      }
    }
  }
  // smoothing
  tmp1(ijk)=F(ijk);
  this->xchng_sclr(tmp1, i, j, k); 
  F(i, j) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j, k) + tmp1(i, j, k - 1));
}

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::surf_sens()
{
  const auto &ijk = this->ijk;
  F(ijk) =  setup::F_sens * hgt_fctr_sclr(ijk); 
// smoothing
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;
  tmp1(ijk)=F(ijk);
  this->xchng_sclr(tmp1, i, j, k); 
  F(i, j, k) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j, k) + tmp1(i, j, k - 1));
}

template <class ct_params_t>
void kin_cloud_3d_lgrngn<ct_params_t>::surf_latent()
{
  const auto &ijk = this->ijk;
  F(ijk) =  setup::F_lat * hgt_fctr_sclr(ijk); 
// smoothing
  const auto &i = this->i;
  const auto &j = this->j;
  const auto &k = this->k;
  tmp1(ijk)=F(ijk);
  this->xchng_sclr(tmp1, i, j, k); 
  F(i, j, k) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j, k) + tmp1(i, j, k - 1));
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
  F(i, j, k) = - w_LS(i, j, k) * (tmp1(i, j, k + 1) - tmp1(i, j, k - 1)) / (2. * this->dk);
  // smoothing
  tmp1(ijk)=F(ijk);
  this->xchng_sclr(tmp1, i, j, k); 
  F(i, j, k) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j, k) + tmp1(i, j, k - 1));
}
