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
  real_t prec_vol;
  std::ofstream f_prec;

  typename parent_t::arr_t w_LS, hgt_fctr_sclr, hgt_fctr_vctr; // TODO: store them in rt_params, here only reference thread's subarrays; also they are just 1D profiles, no need to store whole 3D arrays

  blitz::Array<real_t, 1> k_i; // TODO: make it's size in x direction smaller to match thread's domain

  // global arrays, shared among threads, TODO: in fact no need to share them?
  typename parent_t::arr_t &tmp1,
                           &r_l,
                           &F,       // forcings helper
                           &alpha,   // 'explicit' rhs part - does not depend on the value at n+1
                           &beta;    // 'implicit' rhs part - coefficient of the value at n+1
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

    // recording 1st mom of rd of activated drops
    prtcls->diag_rw_ge_rc();
    prtcls->diag_dry_mom(1);
    this->record_aux("act_rd_mom1", prtcls->outbuf());

    // recording 0th mom of rd of activated drops
    prtcls->diag_rw_ge_rc();
    prtcls->diag_dry_mom(0);
    this->record_aux("act_rd_mom0", prtcls->outbuf());

    // recording 1st mom of rw of gccns
    prtcls->diag_dry_rng(2e-6, 1);
    prtcls->diag_wet_mom(1);
    this->record_aux("gccn_rw_mom1", prtcls->outbuf());

    // recording 0th mom of rw of gccns
    prtcls->diag_dry_rng(2e-6, 1);
    prtcls->diag_wet_mom(0);
    this->record_aux("gccn_rw_mom0", prtcls->outbuf());
   
    // recording total precipitation volume through the lower boundary
    f_prec << this->timestep << " "  << prec_vol << "\n";
    f_prec.flush();
    prec_vol = 0.; 
   
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
    params.cloudph_opts.RH_max = val ? 44 : 1.01; // 0.5% limit during spinup // TODO: specify it somewhere else, dup in blk_2m
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
	  make_arrinfo(this->rhod)
	); 

      // open file for output of precitpitation volume
      f_prec.open(this->outdir+"/prec_vol.dat");
      prec_vol = 0.;

      // writing diagnostic data for the initial condition
      diag();
    }
    // TODO: barrier?
  }

  void vip_rhs_expl_calc()
  {
    parent_t::vip_rhs_expl_calc();
    const auto &i = this->i;
    const auto &j = this->j;
    // kinematic momentum flux  = -u_fric^2 * u_i / |U| * exponential decay
    for (int ji = this->j.first(); ji <= this->j.last(); ++ji)
    {
      F(i, ji) = - pow(setup::u_fric,2) * this->state(ix::vip_i)(i, 0) / sqrt(
                          pow2(this->state(ix::vip_i)(i, 0)))
                          * hgt_fctr_vctr(i, ji);
    }

    // du/dt = sum of kinematic momentum fluxes * dt
    int nz = this->mem->grid_size[1].length(); //76
    blitz::Range notop(0, nz-2);
    this->vip_rhs[0](i, notop) = (F(i, notop) - F(i, notop+1)) / this->dj * this->dt;
    this->vip_rhs[0](i, this->j.last()) = (F(i, this->j.last())) / this->dj * this->dt;
/*
    this->xchng_sclr(F, i, j);
    this->vip_rhs[0](i, j) = (F(i, j) - F(i, j+1)) / this->dj * this->dt;
*/
    // top and bottom cells are two times lower
    this->vip_rhs[0](i, 0) *= 2; 
    this->vip_rhs[0](i, this->j.last()) *= 2; 
  }

  void buoyancy(const blitz::Array<real_t, 2> &th, const blitz::Array<real_t, 2> &rv);
  void radiation(const blitz::Array<real_t, 2> &rv);
  void rv_src();
  void th_src(const blitz::Array<real_t, 2> &rv);
  void w_src(const blitz::Array<real_t, 2> &th, const blitz::Array<real_t, 2> &rv);
  void surf_sens();
  void surf_latent();
  void subsidence(const int&);

  void update_rhs(
    arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  )   
  {   
    parent_t::update_rhs(rhs, dt, at);
    using ix = typename ct_params_t::ix;

    const auto &ijk = this->ijk;

    // forcing
    switch (at) 
    {   
      // for eulerian integration or used to init trapezoidal integration
      case (0): 
      {   
        // ---- water vapor sources ----
        rv_src();
        rhs.at(ix::rv)(ijk) += alpha(ijk) + beta(ijk) * this->state(ix::rv)(ijk); 
        
        // ---- potential temp sources ----
        th_src(this->state(ix::rv));
        rhs.at(ix::th)(ijk) += alpha(ijk) + beta(ijk) * this->state(ix::th)(ijk); 

        // vertical velocity sources
        w_src(this->state(ix::th), this->state(ix::rv));
        rhs.at(ix::w)(ijk) += alpha(ijk);

        // horizontal velocity sources 
        // large-scale vertical wind
        for(auto type : std::set<int>{ix::u})
        {
          subsidence(type);
          rhs.at(type)(ijk) += F(ijk);
        }
        break;
      }   
      case (1): 
      // trapezoidal rhs^n+1
      {   
        // ---- water vapor sources ----
        rv_src();
        rhs.at(ix::rv)(ijk) += (alpha(ijk) + beta(ijk) * this->state(ix::rv)(ijk)) / (1. - 0.5 * this->dt * beta(ijk)); 
        
        // ---- potential temp sources ----
        beta(ijk) = this->state(ix::rv)(ijk) + 0.5 * this->dt * rhs.at(ix::rv)(ijk);
        th_src(beta);
        rhs.at(ix::th)(ijk) += (alpha(ijk) + beta(ijk) * this->state(ix::th)(ijk)) / (1. - 0.5 * this->dt * beta(ijk)); 

        // vertical velocity sources
        // temporarily use beta to store the th^n+1 estimate
        beta(ijk) = this->state(ix::th)(ijk) + 0.5 * this->dt * rhs.at(ix::th)(ijk);
        // temporarily use F to store the rv^n+1 estimate
        F(ijk) = this->state(ix::rv)(ijk) + 0.5 * this->dt * rhs.at(ix::rv)(ijk);
        w_src(beta, F);
        rhs.at(ix::w)(ijk) += alpha(ijk);

        // horizontal velocity sources 
        // large-scale vertical wind
        for(auto type : std::set<int>{ix::u})
        {
          subsidence(type);
          rhs.at(type)(ijk) += F(ijk);
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

    // store -rv_old in tmp1
    tmp1(this->ijk) = - this->mem->advectee(ix::rv)(this->ijk);
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
        prec_vol += ftr.get();
      } else assert(!ftr.valid()); 
#endif

      // store liquid water content to be used in update_rhs (if done in update_rhs, it fails on async runs)
      int nx = this->mem->grid_size[0].length(); //76
      int nz = this->mem->grid_size[1].length(); //76
      prtcls->diag_all();
      prtcls->diag_wet_mom(3);
      auto rl = r_l(blitz::Range(0,nx-1), blitz::Range(0,nz-1)); 
      rl = blitz::Array<real_t,2>(prtcls->outbuf(), blitz::shape(nx, nz), blitz::duplicateData); // copy in data from outbuf; total liquid third moment of wet radius per kg of dry air [m^3 / kg]
      rl = rl * 4./3. * 1000. * 3.14159; // get mixing ratio [kg/kg]
      // in radiation parametrization we integrate mixing ratio * this->rhod
      rl = rl * this->rhod;
      rl = rl * setup::heating_kappa;

      {
        using libmpdataxx::arakawa_c::h;
        // temporarily Cx & Cz are multiplied by this->rhod ...
        auto 
          Cx = this->mem->GC[0](
            this->mem->grid_size[0]^h, 
            this->mem->grid_size[1]
          ).reindex({0,0}).copy(),
          Cz = this->mem->GC[1](
            this->mem->grid_size[0], 
            this->mem->grid_size[1]^h
          ).reindex({0,0}).copy();

        // ... and now dividing them by this->rhod (TODO: z=0 is located at k=1/2)
        {
          blitz::Range all = blitz::Range::all();
          Cx(blitz::Range(1,nx), all)/= this->rhod;
          Cz(all, blitz::Range(1,nz))/= this->rhod;
          Cx(0, all) /= this->rhod(0, all);
          Cz(all, 0) /= this->rhod(all, 0);
        }
        // running synchronous stuff
        prtcls->step_sync(
          params.cloudph_opts,
          make_arrinfo(this->mem->advectee(ix::th)),
          make_arrinfo(this->mem->advectee(ix::rv)),
          make_arrinfo(this->rhod),
          make_arrinfo(Cx), // ix::u ?
          libcloudphxx::lgrngn::arrinfo_t<real_t>(),
          make_arrinfo(Cz) // ix:w ?
        );
        // artificially remove negative rv...
        // this->mem->advectee(ix::rv) = where(this->mem->advectee(ix::rv) < 0., 0., this->mem->advectee(ix::rv));
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
          prec_vol += prtcls->step_async(params.cloudph_opts);
      }
      // performing diagnostics
      if (this->timestep % this->outfreq == 0) 
      { 
#if defined(STD_FUTURE_WORKS)
        if (params.async)
        {
          assert(ftr.valid());
          prec_vol += ftr.get();
        }
#endif
        diag();
      }
    }
    this->mem->barrier();
    // get delta_rv in tmp1
    tmp1(this->ijk) += this->mem->advectee(ix::rv)(this->ijk);
    // update theta accordin to drv
    using libcloudphxx::common::moist_air::c_pd;
    using libcloudphxx::common::const_cp::l_tri;
    // TODO: reenable it after making sure its right
    /*
    for(int x = this->i.first() ; x <= this->i.last(); ++x)
    {
      for(int z = this->j.first() ; z <= this->j.last(); ++z)
      {
        this->mem->advectee(ix::th)(x,z) += - tmp1(x,z) * (l_tri<setup::real_t>() / c_pd<setup::real_t>() / si::kelvin) *
          this->mem->advectee(ix::th)(x,z) / (libcloudphxx::common::theta_dry::T<real_t>(this->state(ix::th)(x, z) * si::kelvins, this->rhod(x, z) * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins); 
      }
    }
*/
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    int backend = -1;
    bool async = true;
    libcloudphxx::lgrngn::opts_t<real_t> cloudph_opts;
    libcloudphxx::lgrngn::opts_init_t<real_t> cloudph_opts_init;
    outmom_t<real_t> out_dry, out_wet;
    real_t z_rlx_sclr;
  };

  private:

  // per-thread copy of params
  rt_params_t params;

  public:

  // ctor
  kin_cloud_2d_lgrngn( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    params(p),
    tmp1(args.mem->tmp[__FILE__][0][0]),
    r_l(args.mem->tmp[__FILE__][0][2]),
    alpha(args.mem->tmp[__FILE__][0][3]),
    beta(args.mem->tmp[__FILE__][0][4]),
    F(args.mem->tmp[__FILE__][0][1])
  {
    int nx = this->mem->grid_size[0].length();
    int nz = this->mem->grid_size[1].length();
    w_LS.resize(nx,nz);
    hgt_fctr_sclr.resize(nx,nz);
    hgt_fctr_vctr.resize(nx,nz);
    k_i.resize(nx);
    r_l = 0.;

    blitz::secondIndex k;

    // prescribed large-scale vertical wind
    w_LS = setup::w_LS_fctr()(k * params.dz);

    // exponential decay with height to distribute constant surface fluxes
    // used to get flux through the bottom of the cell, z=0 at k=1/2
    real_t z_0 = setup::z_rlx_vctr / si::metres;
    hgt_fctr_vctr = exp(- (k-0.5) * params.dz / z_0);
    hgt_fctr_vctr(blitz::Range::all(),0) = 1;
    z_0 = params.z_rlx_sclr;
    hgt_fctr_sclr = exp(- (k-0.5) * params.dz / z_0);
    hgt_fctr_sclr(blitz::Range::all(),0) = 1;

    // delaying any initialisation to ante_loop as rank() does not function within ctor! // TODO: not anymore!!!
    // TODO: equip rank() in libmpdata with an assert() checking if not in serial block
  }  

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 6); // tmp1, tmp2, r_l, alpha, beta, F
  }
};
