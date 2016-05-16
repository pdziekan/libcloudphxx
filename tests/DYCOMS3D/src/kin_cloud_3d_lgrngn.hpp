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
  typename parent_t::arr_t rhod, w_LS, th_init, rv_init, F;
  blitz::Array<real_t, 2> k_i;

  // global arrays, shared among threads
  typename parent_t::arr_t &tmp1,
                           &tmp2,
                           &r_l,
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

  void vip_rhs_expl_calc()
  {
    parent_t::vip_rhs_expl_calc();
    real_t z_0 = setup::z_rlx / si::metres;

    for (int k = this->k.first(); k <= this->k.last(); ++k)
    {
      this->vip_rhs[0](this->i, this->j, k) += - pow(setup::u_fric,2) / z_0 * this->dt / sqrt(
                                                pow2(this->state(ix::vip_i)(this->i, this->j, 0))
                                              + pow2(this->state(ix::vip_j)(this->i, this->j, 0))
                                              ) * this->state(ix::vip_i)(this->i, this->j, 0)
                                                * exp(-this->dk * k / z_0);
      
      this->vip_rhs[1](this->i, this->j, k) += - pow(setup::u_fric,2) / z_0 * this->dt / sqrt(
                                                pow2(this->state(ix::vip_i)(this->i, this->j, 0))
                                              + pow2(this->state(ix::vip_j)(this->i, this->j, 0))
                                              ) * this->state(ix::vip_j)(this->i, this->j, 0)
                                                * exp(-this->dk * k / z_0);
    }
  }

  void buoyancy(const blitz::Array<real_t, 3> &th);
  void radiation(const blitz::Array<real_t, 3> &rv);
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
    const auto &i = this->i;
    const auto &j = this->j;
    const auto &k = this->k;


    // forcing
    switch (at) 
    {   
      // for eulerian integration or used to init trapezoidal integration
      case (0): 
      {   
        // ---- water vapor sources ----
        // surface flux
        surf_latent();
        alpha(i, j, k) = F(i, j, k);
        // large-scale vertical wind
        subsidence(ix::rv);
        alpha(ijk) += tmp1(ijk);
        // TODO: add absorber and nudging to alpha
        //beta(ijk) = 0.;
        // TODO: add absorber and nudging to beta
        rhs.at(ix::rv)(ijk) += alpha(ijk); // TODO: once beta is non-zero, make it alpha + beta * rv
        
        // ---- potential temp sources ----
        // -- heating --
        // surface flux
        surf_sens();
        alpha(i, j, k) = F(i, j, k);
        // radiation
        radiation(this->state(ix::rv));
        alpha(i, j, k) += F(i, j, k);
        // change of theta[K/s] = heating[W/m^3] * theta[K] / T[K] / c_p[J/K/kg] / rhod[kg/m^3]
        for(int x = i.first() ; x <= i.last(); ++x)
        {
          for(int y = j.first() ; y <= j.last(); ++y)
          {
            for(int z = k.first() ; z <= k.last(); ++z)
            {
              alpha(x, y, z) = alpha(x, y, z) * this->state(ix::th)(x, y, z) / rhod(x, y, z) / 
                           (libcloudphxx::common::moist_air::c_p<real_t>(this->state(ix::rv)(x, y, z)) * si::kilograms * si::kelvins / si::joules) / 
                           (libcloudphxx::common::theta_dry::T<real_t>(this->state(ix::th)(x, y, z) * si::kelvins, rhod(x, y, z) * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins);
            }
          }
        }
      
        // large-scale vertical wind
        subsidence(ix::th);
        alpha(ijk) += tmp1(ijk);
        // TODO: add absorber and nudging to alpha
        //beta(ijk) = 0.;
        // TODO: add absorber and nudging to beta
        rhs.at(ix::th)(ijk) += alpha(ijk); // TODO: once beta is non-zero, make it alpha + beta * rv


        // vertical velocity sources
        // buoyancy
        buoyancy(this->state(ix::th));
        alpha(ijk) = tmp2(ijk);
        // large-scale vertical wind
        subsidence(ix::w);
        alpha(ijk) += tmp1(ijk);
        rhs.at(ix::w)(ijk) += alpha(ijk);

        // horizontal velocity sources 
        // large-scale vertical wind
        for(auto type : std::set<int>{ix::u, ix::v})
        {
          subsidence(type);
          rhs.at(type)(ijk) += tmp1(ijk);
        }
        break;
      }   
      case (1): 
      // trapezoidal rhs^n+1
      {   
        // ---- water vapor sources ----
        // surface flux
        surf_latent();
        alpha(i, j, k) = F(i, j, k);
        // large-scale vertical wind
        subsidence(ix::rv);
        alpha(ijk) += tmp1(ijk);
        // TODO: add absorber and nudging to alpha
        //beta(ijk) = 0.;
        // TODO: add absorber and nudging to beta
        rhs.at(ix::rv)(ijk) += alpha(ijk); // TODO: once beta is non-zero, make it (alpha + beta * rv) / (1 - 0.5 * this->dt * beta)
        
        // ---- potential temp sources ----
        // -- heating --
        // surface flux
        surf_sens();
        alpha(i, j, k) = F(i, j, k);
        // temporarily use beta to store the rv^n+1 estimate
        beta(ijk) = this->state(ix::rv)(ijk) + 0.5 * this->dt * rhs.at(ix::rv)(ijk);
        // radiation
        radiation(beta);
        alpha(i, j, k) += F(i, j, k);
        // change of theta[K/s] = heating[W/m^3] * theta[K] / T[K] / c_p[J/K/kg] / rhod[kg/m^3]
        for(int x = i.first() ; x <= i.last(); ++x)
        {
          for(int y = j.first() ; y <= j.last(); ++y)
          {
            for(int z = k.first() ; z <= k.last(); ++z)
            {
              alpha(x, y, z) = alpha(x, y, z) * this->state(ix::th)(x, y, z) / rhod(x, y, z) / 
                           (libcloudphxx::common::moist_air::c_p<real_t>(beta(x, y, z)) * si::kilograms * si::kelvins / si::joules) / 
                           (libcloudphxx::common::theta_dry::T<real_t>(this->state(ix::th)(x, y, z) * si::kelvins, rhod(x, y, z) * si::kilograms / si::metres  / si::metres / si::metres) / si::kelvins);
            }
          }
        }
      
        // large-scale vertical wind
        subsidence(ix::th);
        alpha(ijk) += tmp1(ijk);
        // TODO: add absorber and nudging to alpha
        //beta(ijk) = 0.;
        // TODO: add absorber and nudging to beta
        rhs.at(ix::th)(ijk) += alpha(ijk); // TODO: once beta is non-zero, make it  (alpha + beta * th) / (1 - 0.5 * this->dt * beta)

        // vertical velocity sources
        // temporarily use beta to store the th^n+1 estimate
        beta(ijk) = this->state(ix::th)(ijk) + 0.5 * this->dt * rhs.at(ix::th)(ijk);
        // buoyancy
        buoyancy(beta);
        alpha(ijk) = tmp2(ijk);
        // large-scale vertical wind
        subsidence(ix::w);
        alpha(ijk) += tmp1(ijk);
        rhs.at(ix::w)(ijk) += alpha(ijk);

        // horizontal velocity sources 
        // large-scale vertical wind
        for(auto type : std::set<int>{ix::u, ix::v})
        {
          subsidence(type);
          rhs.at(type)(ijk) += tmp1(ijk);
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

      // store liquid water content to be used in update_rhs (if done in update_rhs, it fails on async runs)
      int nx = this->mem->grid_size[0].length(); //76
      int ny = this->mem->grid_size[1].length(); //76
      int nz = this->mem->grid_size[2].length(); //76
      prtcls->diag_all();
      prtcls->diag_wet_mom(3);
      auto rl = r_l(blitz::Range(0,nx-1), blitz::Range(0,ny-1), blitz::Range(0,nz-1)); 
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
          // rhod is uniformly =1 in mpdata...
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
        // artificially remove negative rv...
        this->mem->advectee(ix::rv) = where(this->mem->advectee(ix::rv) < 0., 0., this->mem->advectee(ix::rv));
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
    r_l(args.mem->tmp[__FILE__][0][2]),
    alpha(args.mem->tmp[__FILE__][0][3]),
    beta(args.mem->tmp[__FILE__][0][4])
  {
    int nx = this->mem->grid_size[0].length();
    int ny = this->mem->grid_size[1].length();
    int nz = this->mem->grid_size[2].length();
    rhod.resize(nx,ny,nz);
    w_LS.resize(nx,ny,nz);
    th_init.resize(nx,ny,nz);
    rv_init.resize(nx,ny,nz);
    F.resize(nx,ny,nz);
    k_i.resize(nx, ny);
    r_l = 0.;

    blitz::thirdIndex k;
    // prescribed density
    rhod = 1;//setup::rhod_fctr()(k * params.dz);

    // prescribed large-scale vertical wind
    w_LS = setup::w_LS_fctr()(k * params.dz);

    // prescribed initial temp profile
    th_init = setup::th_dry_fctr()(k * params.dz);

    // prescribed initial rv profile
    rv_init = 0.; // initially the reference state is not known, will be saved after spinup

    // delaying any initialisation to ante_loop as rank() does not function within ctor! // TODO: not anymore!!!
    // TODO: equip rank() in libmpdata with an assert() checking if not in serial block
  }  

  static void alloc(typename parent_t::mem_t *mem, const int &n_iters)
  {
    parent_t::alloc(mem, n_iters);
    parent_t::alloc_tmp_sclr(mem, __FILE__, 5); // tmp1, tmp2, r_l, alpha, beta
  }
};
