#pragma once

#include "kin_cloud_2d_common.hpp"
#include "outmom.hpp"

#include <libcloudph++/lgrngn/factory.hpp>

#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

// @brief a minimalistic kinematic cloud model with lagrangian microphysics
//        built on top of the mpdata_2d solver (by extending it with
//        custom hook_ante_loop() and hook_post_step() methods)
template <class ct_params_t>
class kin_cloud_2d_lgrngn : public kin_cloud_2d_common<ct_params_t>
{
  // note: lgrngn has no rhs terms - just adjustments (but there might be extrinsic rhs terms)
  using parent_t = kin_cloud_2d_common<ct_params_t>; 

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  private:

  // member fields
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto_t<real_t>> prtcls;

  // helper methods
  void diag()
  {
    assert(this->rank == 0);

    // recording super-droplet concentration per grid cell 
    prtcls->diag_sd_conc();
    this->record_aux("sd_conc", prtcls->outbuf());
   
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
    params.cloudph_opts.coal = params.cloudph_opts.sedi = val; 
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
      if (params.backend != libcloudphxx::lgrngn::CUDA) params.async = false;

      params.cloudph_opts_init.dt = params.dt; // advection timestep = microphysics timestep
      params.cloudph_opts_init.dx = params.dx;
      params.cloudph_opts_init.dz = params.dz;


      // libmpdata++'s grid interpretation
      params.cloudph_opts_init.x0 = params.dx / 2;
      params.cloudph_opts_init.z0 = params.dz / 2;
      params.cloudph_opts_init.x1 = (this->mem->grid_size[0].length() - .5) * params.dx;
      params.cloudph_opts_init.z1 = (this->mem->grid_size[1].length() - .5) * params.dz;

      prtcls.reset(libcloudphxx::lgrngn::factory<real_t>(
        (libcloudphxx::lgrngn::backend_t)params.backend, 
        params.cloudph_opts_init
      ));

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
    parent_t::update_rhs(rhs, dt, at);
    using ix = typename ct_params_t::ix;

    const auto &Tht = this->state(ix::th); 
    const auto &ijk = this->ijk;
    const auto &i = this->i;
    const auto &j = this->j;


    const real_t g = 9.8; 

    switch (at) 
    {   
      case (0): 
      {   
//        rhs.at(ix::tht)(ijk) += H(ijk) - (*this->mem->vab_coeff)(ijk) * (Tht(ijk) - th_init(ijk));

// TODO: update rhs fails for non-serial libmpdata
        tmp1(ijk) = g * (Tht(ijk) - th_init(ijk)) / th_init(ijk);
//std::cout << tmp1 << std::endl;
        this->xchng_sclr(tmp1, i, j); 
//std::cout << tmp1 << std::endl;
        tmp2(i, j) = 0.25 * (tmp1(i, j + 1) + 2 * tmp1(i, j) + tmp1(i, j - 1));
        rhs.at(ix::w)(ijk) += tmp2(ijk);
        // large-scale vertical wind
        for(auto type : std::set<int>{ix::th, ix::rv, ix::u, ix::w})
        {
          tmp1(ijk) = this->state(type)(ijk);
          this->xchng_sclr(tmp1, i, j); 
          tmp2(i, j) = - w_LS(i, j) * (tmp1(i, j + 1) - tmp1(i, j - 1)) / (2. * this->dj);
          rhs.at(type)(ijk) += tmp2(ijk);
        }
        break;
      }   
   
      case (1): 
      {   
//        rhs.at(ix::tht)(ijk) += (H(ijk) - (*this->mem->vab_coeff)(ijk) * (Tht(ijk) - th_init(ijk)))
//                                / (1.0 + 0.5 * (*this->mem->vab_coeff)(ijk) * this->dt);

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
        // --- radiative heating ---
        // TODO: adapt it to trapezoidal integration
        {
          // calc liquid water specific mixing ratio
          int nx = this->mem->grid_size[0].length(); 
          int nz = this->mem->grid_size[1].length(); 
          if(this->rank==0)
          {
//std::cout<<"diag all" << std::endl;
            prtcls->diag_all();
//std::cout<<"diag wet mom 3" << std::endl;
            prtcls->diag_wet_mom(3);
//std::cout<<"diag Arr fron ourbug" << std::endl;
            r_l = blitz::Array<real_t,2>(prtcls->outbuf(), blitz::shape(nx, nz), blitz::neverDeleteData);
//            std::cout << "r_l("<<this->rank<<"): " <<  r_l << std::endl;
            r_l = r_l * 4./3. * 1000. * 3.14159;
            // now multiply by rhod and kappa
            r_l = r_l * rhod; 
            r_l = r_l * setup::heating_kappa;
          }
          this->mem->barrier();
   //       std::cout << "r_l: " <<  r_l << std::endl;
          // find inversion height
          blitz::secondIndex j;
          j_i = blitz::first( (this->mem->advectee(ix::rv) + r_l) < setup::q_i, j);
     //     std::cout << "j_i: " << j_i << std::endl;

          // calc Eqs. 5 and 6 from Ackerman et al 2009
          for(int i = 0 ; i < nx; ++i)
          {
            for(int j = 0 ; j < nz; ++j)
            {
              Q(i, j) = this->mem->sum(r_l, rng_t(i, i), rng_t(j, this->j.last()), false);
              Q(i, j) = Q(i, j) * (nz - j) * this->dj;
              F(i, j) = setup::F_0 * exp(- Q(i, j));

              if(j > 0)
              {
                Q(i, j) = this->mem->sum(r_l, rng_t(i, i), rng_t(0, j-1), false);
                Q(i, j) = Q(i, j) * (nz - j) * this->dj;
                F(i, j) += setup::F_0 * exp(- Q(i, j));
              }
              if(j > j_i(i) )
              {
                real_t z_d = (j - j_i(i)) * this->dj;
                F(i, j) += setup::c_p * setup::rho_i * setup::D * (0.25 * pow(z_d, 4./3) + j_i(i) * this->dj * pow(z_d, 1./3)); 
              }
     //         std::cout << "F: " << F << std::endl;
              F(i, j) = F(i, j) / (libcloudphxx::common::moist_air::c_p<real_t>(this->state(ix::rv)(i, j)) * si::kilograms * si::kelvins / si::joules); // divide by specific heat capacity
            }
          }
          rhs.at(ix::th)(ijk) += F(ijk) / this->dj / rhod(ijk); // heating[W/m^2] / cell height[m] / rhod[kg/m^3] / specific heat capacity of moist air [J/K/kg]
        }
        // --- surface fluxes ---
        {
          int nx = this->mem->grid_size[0].length(); 

          // sensible heat
          for(int i = 0 ; i < nx; ++i)
          {
            F(i, 0) = setup::F_sens / (libcloudphxx::common::moist_air::c_p<real_t>(this->state(ix::rv)(i, 0)) * si::kilograms * si::kelvins / si::joules); // divide by specific heat capacity
          }
          rhs.at(ix::th)(i, blitz::Range(0,0)) += F(i, blitz::Range(0,0)) / this->dj / rhod(i, blitz::Range(0,0)); // heating[W/m^2] / cell height[m] / rhod[kg/m^3] / specific heat capacity of moist air [J/K/kg]

          // latent heat
          rhs.at(ix::rv)(i, blitz::Range(0,0)) += setup::F_lat / (libcloudphxx::common::const_cp::l_tri<real_t>() * si::kilograms / si::joules) / this->dj / rhod(i, blitz::Range(0,0)); // heating[W/m^2] / cell height[m] / rhod[kg/m^3] / latent heat of evaporation [J/kg]
 
          // momentum flux
          for(int i = 0 ; i < nx; ++i)
          {
            F(i, 0) = this->state(ix::u)(i, 0) < 0 ? 1 : -1;  // sign of u
          }
          rhs.at(ix::u)(i, blitz::Range(0,0)) += F(i, blitz::Range(0,0)) *  pow(setup::u_fric,2) /  this->dj; 
        }
        break;
      }   
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

      {
        using libmpdataxx::arakawa_c::h;
        // temporarily Cx & Cz are multiplied by rhod ...
        auto 
          Cx = this->mem->GC[0](
            this->mem->grid_size[0]^h, 
            this->mem->grid_size[1]
          ).reindex({0,0}).copy(),
          Cz = this->mem->GC[1](
            this->mem->grid_size[0], 
            this->mem->grid_size[1]^h
          ).reindex({0,0}).copy();

        // ... and now dividing them by rhod (z=0 is located at j=1/2)
        {
          blitz::secondIndex j;
          Cx /= setup::rhod()(   j     * this->dj);
          Cz /= setup::rhod()((j - .5) * this->dj);
        }
        // running synchronous stuff
        prtcls->step_sync(
          params.cloudph_opts,
          make_arrinfo(this->mem->advectee(ix::th)),
          make_arrinfo(this->mem->advectee(ix::rv)),
          make_arrinfo(rhod),
          make_arrinfo(Cx), // ix::u ?
          libcloudphxx::lgrngn::arrinfo_t<real_t>(),
          make_arrinfo(Cz) // ix:w ?
        );
      } 

      // running asynchronous stuff
      {
        using libcloudphxx::lgrngn::particles_t;
        using libcloudphxx::lgrngn::CUDA;

#if defined(STD_FUTURE_WORKS)
        if (params.async)
        {
          assert(!ftr.valid());
          ftr = std::async(
            std::launch::async, 
            &particles_t<real_t, CUDA>::step_async, 
            dynamic_cast<particles_t<real_t, CUDA>*>(prtcls.get()),
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
  blitz::Array<real_t,2> th_init;
  blitz::Array<real_t,2> tmp1;
  blitz::Array<real_t,2> tmp2;
  blitz::Array<real_t,2> rhod;
  blitz::Array<real_t,2> w_LS; //large-scale vertical wind
  blitz::Array<real_t,2> F; // radiatvie heating
  blitz::Array<real_t,2> Q; // radiatvie heating hlpr
  blitz::Array<real_t,1> j_i; // inversion height index
  blitz::Array<real_t,2> r_l; 

  public:

  // ctor
  kin_cloud_2d_lgrngn( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    params(p)
  {
    th_init.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    rhod.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    tmp1.resize(blitz::Range(-2, this->mem->grid_size[0].length() + 2), blitz::Range(-2, this->mem->grid_size[1].length() + 2));
    tmp2.resize(blitz::Range(-2, this->mem->grid_size[0].length() + 2), blitz::Range(-2, this->mem->grid_size[1].length() + 2));
    w_LS.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    F.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    Q.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    r_l.resize(this->mem->grid_size[0].length(), this->mem->grid_size[1].length());
    j_i.resize(this->mem->grid_size[0].length());
    blitz::secondIndex j;
    th_init = setup::th_dry()(j * p.dz);
    rhod = setup::rhod()(j * p.dz);
    w_LS = setup::w_LS()(j * p.dz);

    // delaying any initialisation to ante_loop as rank() does not function within ctor! // TODO: not anymore!!!
    // TODO: equip rank() in libmpdata with an assert() checking if not in serial block
  }  
};
