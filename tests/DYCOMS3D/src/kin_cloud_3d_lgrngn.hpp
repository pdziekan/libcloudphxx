#pragma once

#include "kin_cloud_3d_common.hpp"
#include "outmom.hpp"

#include <libcloudph++/lgrngn/factory.hpp>

#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

// @brief a minimalistic kinematic cloud model with lagrangian microphysics
//        built on top of the mpdata_2d solver (by extending it with
//        custom hook_ante_loop() and hook_post_step() methods)
template <class ct_params_t>
class kin_cloud_3d_lgrngn : public kin_cloud_3d_common<ct_params_t>
{
  // note: lgrngn has no rhs terms - just adjustments (but there might be extrinsic rhs terms)
  using parent_t = kin_cloud_3d_common<ct_params_t>; 

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  private:

  // member fields
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto_t<real_t>> prtcls;
  typename parent_t::arr_t rhod, w_LS, th_init, F;
  blitz::Array<real_t, 2> k_i;

  // global arrays, shared among threads
  typename parent_t::arr_t &tmp1,
                           &tmp2,
                           &r_l;  // these 3 are modified during simulation
/*
                           &rhod,
                           &w_LS, 
                           &th_init, // these 3 are set during init
*/
/*                           &F, // radiatvie heating
                           &Q, // radiatvie heating hlpr
                           &j_i; // inversion height index */
  // helper methods
  void diag()
  {
    assert(this->rank == 0);

    // recording super-droplet concentration per grid cell 
    prtcls->diag_sd_conc();
    this->record_aux("sd_conc", prtcls->outbuf());

    // recording precipitation rate per grid cel
    prtcls->diag_all();
    prtcls->diag_precip_rate();
    this->record_aux("precip_rate", prtcls->outbuf());
   
    // recording requested statistical moments
    {
      // dry
      int rng_num = 0;
      for (auto &rng_moms : params.out_dry)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_dry_rng(rng.first / si::metres, rng.second / si::metres);
        for (auto &mom : rng_moms.second)
        {
          prtcls->diag_dry_mom(mom);
          this->record_aux(aux_name("rd", rng_num, mom), prtcls->outbuf());
        }
        rng_num++;
      }
    }
    {
      // wet
      int rng_num = 0;
      for (auto &rng_moms : params.out_wet)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_wet_rng(rng.first / si::metres, rng.second / si::metres);
        for (auto &mom : rng_moms.second)
        {
          prtcls->diag_wet_mom(mom);
          this->record_aux(aux_name("rw", rng_num, mom), prtcls->outbuf());
        }
        rng_num++;
      }
    }
  } 

  libcloudphxx::lgrngn::arrinfo_t<real_t> make_arrinfo(
    typename parent_t::arr_t arr
  ) {
    return libcloudphxx::lgrngn::arrinfo_t<real_t>(
      arr.dataZero(), 
      arr.stride().data()
    );
  }

  std::string aux_name(
    const std::string pfx, 
    const int rng,
    const int mom
  )
  { 
    std::ostringstream tmp;
    tmp << pfx << "_rng" << std::setw(3) << std::setfill('0') << rng << "_mom" << mom;
    return tmp.str();
  }

  protected:

  bool get_rain() { return params.cloudph_opts.coal && params.cloudph_opts.sedi; }
  void set_rain(bool val) 
  { 
    params.cloudph_opts.coal = val;
    params.cloudph_opts.sedi = val; 
    params.cloudph_opts.RH_max = val ? 44 : 1.01; // 1% limit during spinup // TODO: specify it somewhere else, dup in blk_2m
  };

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    parent_t::hook_ante_loop(nt); 

    // TODO: barrier?
    if (this->rank == 0) 
    {
      assert(params.backend != -1);
      assert(params.dt != 0); 

      // async does not make sense without CUDA
      if (params.backend != libcloudphxx::lgrngn::CUDA && params.backend != libcloudphxx::lgrngn::multi_CUDA) params.async = false;

      params.cloudph_opts_init.dt = params.dt; // advection timestep = microphysics timestep
      params.cloudph_opts_init.dx = params.dx;
      params.cloudph_opts_init.dy = params.dy;
      params.cloudph_opts_init.dz = params.dz;


      // libmpdata++'s grid interpretation
      params.cloudph_opts_init.x0 = params.dx / 2;
      params.cloudph_opts_init.y0 = params.dy / 2;
      params.cloudph_opts_init.z0 = params.dz / 2;
      params.cloudph_opts_init.x1 = (this->mem->grid_size[0].length() - .5) * params.dx;
      params.cloudph_opts_init.y1 = (this->mem->grid_size[1].length() - .5) * params.dy;
      params.cloudph_opts_init.z1 = (this->mem->grid_size[2].length() - .5) * params.dz;

      prtcls.reset(libcloudphxx::lgrngn::factory<real_t>(
        (libcloudphxx::lgrngn::backend_t)params.backend, 
        params.cloudph_opts_init
      ));

//  std::cout << "call init" << std::endl;
//  std::cout << "rhod " << rhod << std::endl;
	prtcls->init(
	  make_arrinfo(this->mem->advectee(ix::th)),
	  make_arrinfo(this->mem->advectee(ix::rv)),
	  make_arrinfo(rhod)
	); 
      // writing diagnostic data for the initial condition
      diag();
    }
    // TODO: barrier?
  }

  void update_rhs(
    arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  )   
  {   
//  std::cout << "call updatre rhs" << std::endl;
    parent_t::update_rhs(rhs, dt, at);
    using ix = typename ct_params_t::ix;

    const auto &Tht = this->state(ix::th); 
    const auto &rv = this->state(ix::rv); 
    const auto &ijk = this->ijk;
    const auto &i = this->i;
    const auto &j = this->j;
    const auto &k = this->k;

    const real_t g = 9.8; 

    // forcing
    switch (at) 
    {   
      case (0): 
      {   
        int nx = this->mem->grid_size[0].length(); //76
        int ny = this->mem->grid_size[1].length(); //76
        int nz = this->mem->grid_size[2].length(); //76
// TODO: update rhs fails for non-serial libmpdata
        // buoyancy
        tmp1(ijk) = g * (Tht(ijk) - th_init(ijk)) / th_init(ijk);
        this->xchng_sclr(tmp1, i, j, k); 
        tmp2(i, j, k) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j, k) + tmp1(i, j, k - 1));
        rhs.at(ix::w)(ijk) += tmp2(ijk);

        // large-scale vertical wind
        for(auto type : std::set<int>{ix::th, ix::rv, ix::u, ix::w})
        {
          tmp1(ijk) = this->state(type)(ijk);
          this->xchng_sclr(tmp1, i, j, k);
          tmp2(i, j, k) = - w_LS(i, j, k) * (tmp1(i, j, k + 1) - tmp1(i, j, k - 1)) / (2. * this->dk); // use built-in blitz stencil?
          tmp1(i, j, k) = 0.25 * (tmp2(i, j, k + 1) + 2 * tmp2(i, j, k) + tmp2(i, j, k - 1));
          rhs.at(type)(ijk) += tmp1(ijk);
        }
        // --- radiative heating ---
        // TODO: adapt it to trapezoidal integration
        {
//          std::cout << "nx: " << nx << "nz: " << nz << std::endl;

          // index of first cell above inversion
          blitz::thirdIndex ki;
          blitz::Array<real_t, 3> tmp(rv(ijk).shape());
          tmp  = rv(ijk) + r_l(ijk);
          k_i(i, j) = blitz::first( tmp< setup::q_i, ki) - rv.base(blitz::secondDim); // rv and r_l have same bases (first indices), subarrays (i.e. rv(ijk)) start with the same base as original arr (i.e. rv)!

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

                F(x, y, z) = F(x, y, z) / (libcloudphxx::common::moist_air::c_p<real_t>(rv(x, y, z)) * si::kilograms * si::kelvins / si::joules); // divide by specific heat capacity
                F(x, y, z) = F(x, y, z) / (libcloudphxx::common::theta_dry::T<real_t>(Tht(x, y, z) * si::kelvins, rhod(x, y, z) * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins); // divide by temperature
              }
            }
          }
          F(ijk) = F(ijk) / this->dk / rhod(ijk); // heating[W/m^2] / cell height[m] / rhod[kg/m^3] / specific heat capacity of moist air [J/K/kg]
          // Eq. 3.33 from Curry and Webster
          blitz::Range notopbot(1, nz-2);
          tmp1(i, j, notopbot) -= Tht(i, j, notopbot) *                          // theta dry
                             (F(i, j, notopbot+1) - F(i, j, notopbot-1)) / 2.;   // gradient of heat flux 
          tmp1(i, j, k.last()) -= Tht(i, j, k.last()) *                          // theta dry
                               (F(i, j, k.last()) - F(i, j, k.last()-1));       // gradient of heat flux 
          tmp1(i, j, 0) -= Tht(i, j, 0) *                                        // theta dry
                        (F(i, j, 1) - F(i, j, 0));                              // gradient of heat flux 
          this->xchng_sclr(tmp1, i, j, k);
          tmp2(i, j, k) = 0.25 * (tmp1(i, j, k + 1) + 2 * tmp1(i, j, k) + tmp1(i, j, k - 1));
          rhs.at(ix::th)(ijk) += tmp2(ijk);
// debug output
/*
for(int rank=0;rank<4;++rank)
{
if(this->rank==rank)
{
          std::cout << "rv " << rv << std::endl;
          std::cout << "rv(ijk) " << rv(ijk) << std::endl;
          std::cout << "rl " << r_l << std::endl;
          std::cout << "rl(ijk) " << r_l(ijk) << std::endl;
          std::cout << "rhod(ijk) " << rhod(ijk) << std::endl;
          std::cout << "F(ijk) " << F(ijk) << std::endl;
          std::cout << "rhs(th) " << rhs.at(ix::th) << std::endl;
          std::cout << "th " << Tht << std::endl;
//          std::cout << "first " << blitz::first( (rv + r_l ) < setup::q_i, ji) << std::endl;
          std::cout << "ji " << j_i(i) << std::endl;
          std::cout << this->rank << " " << i.first() << " " << i.last() << std::endl;
          std::cout << this->rank << " " << j.first() << " " << j.last() << std::endl;
          std::cout << this->rank << " " << ijk.lbound(0) << " " << ijk.ubound(0) << std::endl;
          std::cout << this->rank << " " << ijk.lbound(1) << " " << ijk.ubound(1) << std::endl;
}
this->mem->barrier();
}
*/
        }
        // --- surface fluxes ---
        {
          // sensible heat
          for(int x = i.first() ; x <= i.last(); ++x)
          {
            for(int y = j.first() ; y <= j.last(); ++y)
            {
              F(x, y, 0) = setup::F_sens / (libcloudphxx::common::moist_air::c_p<real_t>(this->state(ix::rv)(x, y, 0)) * si::kilograms * si::kelvins / si::joules); // heating divided by specific heat capacity
              F(x, y, 0) = F(x, y, 0) / (libcloudphxx::common::theta_dry::T<real_t>(Tht(x, y, 0) * si::kelvins, rhod(x, y, 0) * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins); // divide by temperature
            }
          }
          // heating[W/m^2] / cell height[m] / rhod[kg/m^3] / specific heat capacity of moist air [J/K/kg]
          rhs.at(ix::th)(i, j, blitz::Range(0,0)) += F(i, j, blitz::Range(0,0)) *                                                       // heat in W/m^2 divided by spec heat capacity
            Tht(i, j, blitz::Range(0,0)) /                                                                                           // times dry potential temp
            this->dk /                                                                                                            // divide by cell height
            rhod(i, j, blitz::Range(0,0));                                                                                           // divide by density

          // latent heat
          // heating[W/m^2] / cell height[m] / rhod[kg/m^3] / latent heat of evaporation [J/kg]
          rhs.at(ix::rv)(i, j, blitz::Range(0,0)) += setup::F_lat /                           // heating 
          (libcloudphxx::common::const_cp::l_tri<real_t>() * si::kilograms / si::joules) / // latent heat of evaporation
          this->dk /                                                                       // cell height
          rhod(i, j, blitz::Range(0,0));                                                      // density
 
          // momentum flux
          blitz::Array<real_t, 3> uMag(nx, ny);
          uMag = sqrt(
                   this->state(ix::u)(i, j, blitz::Range(0,0)) * this->state(ix::u)(i, j, blitz::Range(0,0)) +
                   this->state(ix::v)(i, j, blitz::Range(0,0)) * this->state(ix::v)(i, j, blitz::Range(0,0))
                 );
          rhs.at(ix::u)(i, j, blitz::Range(0,0)) -= this->state(ix::u)(i, j, blitz::Range(0,0)) / uMag(i, j) *  pow(setup::u_fric,2) /  this->dk;  
          rhs.at(ix::v)(i, j, blitz::Range(0,0)) -= this->state(ix::v)(i, j, blitz::Range(0,0)) / uMag(i, j) *  pow(setup::u_fric,2) /  this->dk;  
        }
        break;
      }   
/*   
      case (1): 
      {   

        // buoyancy
        tmp1(ijk) = g * (Tht(ijk) + 0.5 * this->dt * rhs.at(ix::th)(ijk) - th_init(ijk)) / th_init(ijk);
        this->xchng_sclr(tmp1, i, j); 
        tmp2(i, j) = 0.25 * (tmp1(i, j + 1) + 2 * tmp1(i, j) + tmp1(i, j - 1));
        rhs.at(ix::w)(ijk) += tmp2(ijk);

        // large-scale vertical wind
        for(auto type : std::set<int>{ix::th, ix::rv, ix::u, ix::w})
        {
          tmp1(ijk) = this->state(type)(ijk) + 0.5 * this->dt * rhs.at(type)(ijk);
          this->xchng_sclr(tmp1, i, j); 
          tmp2(i, j) = - w_LS(i, j) * (tmp1(i, j + 1) - tmp1(i, j - 1)) / (2. * this->dj);
          rhs.at(type)(ijk) += tmp2(ijk);
        }
        break;
      }
*/
    }  
  }


#if defined(STD_FUTURE_WORKS)
  std::future<real_t> ftr;
#endif

  // 
  void hook_post_step()
  {
    parent_t::hook_post_step(); // includes output

    this->mem->barrier();
    if (this->rank == 0) 
    {
      // assuring previous async step finished ...
#if defined(STD_FUTURE_WORKS)
      if (
        params.async && 
        this->timestep != 0 && // ... but not in first timestep ...
        ((this->timestep - 1) % this->outfreq != 0) // ... and not after diag call
      ) {
        assert(ftr.valid());
        ftr.get();
      } else assert(!ftr.valid()); 
#endif

      // store liquid water content to be used in update_rhs (if done in update_rhs, it fails on async runs)
      int nx = this->mem->grid_size[0].length(); //76
      int ny = this->mem->grid_size[1].length(); //76
      int nz = this->mem->grid_size[2].length(); //76
      prtcls->diag_all();
      prtcls->diag_wet_mom(3);
      auto rl = r_l(blitz::Range(0,nx-1), blitz::Range(0,ny-1), blitz::Range(0,nz-1)); // rl references the nohalo subdomain
      rl = blitz::Array<real_t,3>(prtcls->outbuf(), blitz::shape(nx, ny, nz), blitz::duplicateData); // copy in data from outbuf
      rl = rl * 4./3. * 1000. * 3.14159;
      rl = rl * rhod; 
      rl = rl * setup::heating_kappa;

      {
        using libmpdataxx::arakawa_c::h;
        // temporarily Cx & Cz are multiplied by rhod ...
        auto 
          Cx = this->mem->GC[0](
            this->mem->grid_size[0]^h, 
            this->mem->grid_size[1],
            this->mem->grid_size[2]
          ).reindex({0,0,0}).copy(),
          Cy = this->mem->GC[1](
            this->mem->grid_size[0], 
            this->mem->grid_size[1]^h,
            this->mem->grid_size[2]
          ).reindex({0,0,0}).copy(),
          Cz = this->mem->GC[2](
            this->mem->grid_size[0], 
            this->mem->grid_size[1], 
            this->mem->grid_size[2]^h
          ).reindex({0,0,0}).copy();

        // ... and now dividing them by rhod (z=0 is located at k=1/2)
        {
          blitz::thirdIndex k;
          Cx /= setup::rhod_fctr()(   k     * this->dk);
          Cy /= setup::rhod_fctr()(   k     * this->dk);
          Cz /= setup::rhod_fctr()((k - .5) * this->dk);
        }
        // running synchronous stuff
        prtcls->step_sync(
          params.cloudph_opts,
          make_arrinfo(this->mem->advectee(ix::th)),
          make_arrinfo(this->mem->advectee(ix::rv)),
          make_arrinfo(rhod),
          make_arrinfo(Cx), // ix::u ?
          make_arrinfo(Cy), // ix::u ?
          make_arrinfo(Cz) // ix:w ?
        );
      } 

      // running asynchronous stuff
      {
        using libcloudphxx::lgrngn::particles_t;
        using libcloudphxx::lgrngn::CUDA;
        using libcloudphxx::lgrngn::multi_CUDA;

#if defined(STD_FUTURE_WORKS)
        if (params.async)
        {
          assert(!ftr.valid());
          if(params.backend == CUDA)
            ftr = std::async(
              std::launch::async, 
              &particles_t<real_t, CUDA>::step_async, 
              dynamic_cast<particles_t<real_t, CUDA>*>(prtcls.get()),
              params.cloudph_opts
            );
          else if(params.backend == multi_CUDA)
            ftr = std::async(
              std::launch::async, 
              &particles_t<real_t, multi_CUDA>::step_async, 
              dynamic_cast<particles_t<real_t, multi_CUDA>*>(prtcls.get()),
              params.cloudph_opts
            );
          assert(ftr.valid());
        } else 
#endif
          prtcls->step_async(params.cloudph_opts);
      }

      // performing diagnostics
      if (this->timestep % this->outfreq == 0) 
      { 
#if defined(STD_FUTURE_WORKS)
        if (params.async)
        {
          assert(ftr.valid());
          ftr.get();
        }
#endif
        diag();
      }
    }

    this->mem->barrier();
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    int backend = -1;
    bool async = true;
    libcloudphxx::lgrngn::opts_t<real_t> cloudph_opts;
    libcloudphxx::lgrngn::opts_init_t<real_t> cloudph_opts_init;
    outmom_t<real_t> out_dry, out_wet;
  };

  private:

  // per-thread copy of params
  rt_params_t params;

  public:

  // ctor
  kin_cloud_3d_lgrngn( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    params(p),
    tmp1(args.mem->tmp[__FILE__][0][0]),
    tmp2(args.mem->tmp[__FILE__][0][1]),
    r_l(args.mem->tmp[__FILE__][0][2])
  {
    int nx = this->mem->grid_size[0].length();
    int ny = this->mem->grid_size[1].length();
    int nz = this->mem->grid_size[2].length();
    rhod.resize(nx,ny,nz);
    w_LS.resize(nx,ny,nz);
    th_init.resize(nx,ny,nz);
    F.resize(nx,ny,nz);
    k_i.resize(nx, ny);

    blitz::thirdIndex k;
    // prescribed density
    rhod = setup::rhod_fctr()(k * params.dz);
//    std::cout << "rhod w setopts" << rhod << std::endl;

    // prescribed large-scale vertical wind
    w_LS = setup::w_LS_fctr()(k * params.dz);

    // prescribed initial temp profile
    th_init = setup::th_dry_fctr()(k * params.dz);

/*    th_init.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    rhod.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    tmp1.resize(blitz::Range(-2, this->mem->grid_size[0].length() + 2), blitz::Range(-2, this->mem->grid_size[1].length() + 2));
    tmp2.resize(blitz::Range(-2, this->mem->grid_size[0].length() + 2), blitz::Range(-2, this->mem->grid_size[1].length() + 2));
    w_LS.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    F.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    Q.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    r_l.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    j_i.resize(this->mem->grid_size[0].length());
*/

    // delaying any initialisation to ante_loop as rank() does not function within ctor! // TODO: not anymore!!!
    // TODO: equip rank() in libmpdata with an assert() checking if not in serial block
  }  

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 3); // tmp1, tmp2, r_l
  }
};
