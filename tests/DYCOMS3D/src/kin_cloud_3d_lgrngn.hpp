#pragma once

#include "kin_cloud_3d_common.hpp"
#include "outmom.hpp"

#include <libcloudph++/lgrngn/factory.hpp>

#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

#include <chrono>


// @brief a minimalistic kinematic cloud model with lagrangian microphysics
//        built on top of the mpdata_2d solver (by extending it with
//        custom hook_ante_loop() and hook_post_step() methods)
template <class ct_params_t>
class kin_cloud_3d_lgrngn : public kin_cloud_3d_common<ct_params_t>
{
  // note: lgrngn has no rhs terms - just adjustments (but there might be extrinsic rhs terms)
  using parent_t = kin_cloud_3d_common<ct_params_t>; 
  using clock = std::chrono::high_resolution_clock;

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  private:

  // member fields
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto_t<real_t>> prtcls;
  real_t prec_vol;
  std::ofstream f_prec;
  clock::time_point tbeg, tend, tbeg1, tend1, tbeg_loop;
  std::chrono::milliseconds tdiag, tupdate, tsync, tasync_wait, tloop, tpost_step_custom, tpost_step_base;

  typename parent_t::arr_t rhod, w_LS, th_init, rv_init, hgt_fctr_sclr, hgt_fctr_vctr; // TODO: store them in rt_params, here only reference thread's subarrays; also they are just 1D profiles, no need to store whole 3D arrays
  blitz::Array<real_t, 2> k_i; // TODO: make it's size in x direction smaller to match thread's domain

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
    tbeg = clock::now();

    // recording super-droplet concentration per grid cell 
    prtcls->diag_sd_conc();
    this->record_aux("sd_conc", prtcls->outbuf());

    // recording precipitation rate per grid cel
    prtcls->diag_all();
    prtcls->diag_precip_rate();
    this->record_aux("precip_rate", prtcls->outbuf());
   
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
    tend = clock::now();
    tdiag += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
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
//      if (params.backend != libcloudphxx::lgrngn::CUDA && params.backend != libcloudphxx::lgrngn::multi_CUDA) params.async = false;

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
	  make_arrinfo(this->rhod)
	); 

      // open file for output of precitpitation volume
      f_prec.open(this->outdir+"/prec_vol.dat");
      prec_vol = 0.;

      // writing diagnostic data for the initial condition
      diag();
      tbeg_loop = clock::now();
    }
    // TODO: barrier?
  }

  void vip_rhs_expl_calc()
  {
    parent_t::vip_rhs_expl_calc();
    const auto &i = this->i;
    const auto &j = this->j;
    const auto &k = this->k;
    int nz = this->mem->grid_size[2].length(); //76
    for(int d=0; d<2; ++d)
    {
      auto vip_ij = d==0? ix::vip_i : ix::vip_j;
      for (int ki = this->k.first(); ki <= this->k.last(); ++ki)
      {
        // kinematic momentum flux = -u_fric^2 * u_i / |U| * exponential decay
        F(i, j, ki) = - pow(setup::u_fric,2) * this->state(vip_ij)(i, j, 0) / sqrt(
                              pow2(this->state(ix::vip_i)(i, j, 0))) 
                              * hgt_fctr_vctr(i, j, ki);
      }
      // du/dt = sum of kinematic momentum fluxes * dt
      blitz::Range notop(0, nz-2);
      this->vip_rhs[d](i, j, notop) = ( F(i, j, notop) - F(i, j, notop+1)) / this->dk * this->dt;
      this->vip_rhs[d](i, j, k.last()) = ( F(i, j, k.last())) / this->dk * this->dt;

      // top and bottom cells are two times lower
      this->vip_rhs[d](i, j, 0) *= 2; 
      this->vip_rhs[d](i, j, k.last()) *= 2; 
    }
  }

  void buoyancy(const blitz::Array<real_t, 3> &th, const blitz::Array<real_t, 3> &rv);
  void radiation(const blitz::Array<real_t, 3> &rv);
  void rv_src();
  void th_src(const blitz::Array<real_t, 3> &rv);
  void w_src(const blitz::Array<real_t, 3> &th, const blitz::Array<real_t, 3> &rv);
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
    this->mem->barrier();
    if(this->rank == 0)
      tbeg = clock::now();

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
        for(auto type : std::set<int>{ix::u, ix::v})
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
        // temporarily use beta to store the rv^n+1 estimate
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
        for(auto type : std::set<int>{ix::u, ix::v})
        {
          subsidence(type);
          rhs.at(type)(ijk) += F(ijk);
        }
        break;
      }
    }  
    this->mem->barrier();
    if(this->rank == 0)
    {
      tend = clock::now();
      tupdate += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
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
      if (blitz::min(this->mem->advectee(ix::rv)) < 0.)
      {
        std::cout << "timestep " << this->timestep << " negative rv in hook_post_step: " << this->mem->advectee(ix::rv) << std::endl;
        // artificially remove negative rv...
        this->mem->advectee(ix::rv) = where(this->mem->advectee(ix::rv) < 0., 0., this->mem->advectee(ix::rv));
      }
      tbeg1 = clock::now();
      // assuring previous async step finished ...
#if defined(STD_FUTURE_WORKS)
      if (
        params.async && 
        this->timestep != 0 && // ... but not in first timestep ...
        ((this->timestep - 1) % this->outfreq != 0) // ... and not after diag call
      ) {
        assert(ftr.valid());
        tbeg = clock::now();
        prec_vol += ftr.get();
        tend = clock::now();
        tasync_wait += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
      } else assert(!ftr.valid()); 
#endif

      // store liquid water content to be used in update_rhs (if done in update_rhs, it fails on async runs)
      int nx = this->mem->grid_size[0].length(); //76
      int ny = this->mem->grid_size[1].length(); //76
      int nz = this->mem->grid_size[2].length(); //76
      prtcls->diag_all();
      prtcls->diag_wet_mom(3);
      auto rl = r_l(blitz::Range(0,nx-1), blitz::Range(0,ny-1), blitz::Range(0,nz-1)); 
      rl = blitz::Array<real_t,3>(prtcls->outbuf(), blitz::shape(nx, ny, nz), blitz::neverDeleteData); // copy in data from outbuf
      rl = rl * 4./3. * 1000. * 3.14159 * this->rhod * setup::heating_kappa;
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

        // ... and now dividing them by rhod (TODO: z=0 is located at k=1/2)
        {
          blitz::Range all = blitz::Range::all();
          Cx(blitz::Range(1,nx), all, all)/= this->rhod;
          Cy(all, blitz::Range(1,ny), all)/= this->rhod;
          Cz(all, all, blitz::Range(1,nz))/= this->rhod;
          Cx(0, all, all) /= this->rhod(0, all, all);
          Cy(all, 0, all) /= this->rhod(all, 0, all);
          Cz(all, all, 0) /= this->rhod(all, all, 0);
        }
        // running synchronous stuff
        tbeg = clock::now();
        prtcls->step_sync(
          params.cloudph_opts,
          make_arrinfo(this->mem->advectee(ix::th)),
          make_arrinfo(this->mem->advectee(ix::rv)),
          make_arrinfo(this->rhod),
          make_arrinfo(Cx), // ix::u ?
          make_arrinfo(Cy), // ix::u ?
          make_arrinfo(Cz) // ix:w ?
        );
        tend = clock::now();
        tsync += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
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
      tend1 = clock::now();
      tpost_step_custom += std::chrono::duration_cast<std::chrono::milliseconds>( tend1 - tbeg1 );
      // there's no hook_post_loop, so we imitate it here
      if(this->timestep == params.nt-1)
      {
        tend = clock::now();
        tloop = std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg_loop );
        std::cout <<  "wall time in milliseconds: " << std::endl
          << "loop: " << tloop.count() << std::endl
          << "update: " << tupdate.count() << " ("<< setup::real_t(tupdate.count())/tloop.count()*100 <<"%)" << std::endl
          << "custom_post_step: " << tpost_step_custom.count() << " ("<< setup::real_t(tpost_step_custom.count())/tloop.count()*100 <<"%)" << std::endl
          << "base_post_step: " << tpost_step_base.count() << " ("<< setup::real_t(tpost_step_base.count())/tloop.count()*100 <<"%)" << std::endl
          << "diag: " << tdiag.count() << " ("<< setup::real_t(tdiag.count())/tloop.count()*100 <<"%)" << std::endl
          << "sync: " << tsync.count() << " ("<< setup::real_t(tsync.count())/tloop.count()*100 <<"%)" << std::endl
          << "async_wait: " << tasync_wait.count() << " ("<< setup::real_t(tasync_wait.count())/tloop.count()*100 <<"%)" << std::endl;
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
    real_t z_rlx_sclr;
    int nt; // number of timesteps
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
    r_l(args.mem->tmp[__FILE__][0][2]),
    alpha(args.mem->tmp[__FILE__][0][3]),
    beta(args.mem->tmp[__FILE__][0][4]),
    F(args.mem->tmp[__FILE__][0][1]),
    tdiag(std::chrono::milliseconds::zero()),
    tupdate(std::chrono::milliseconds::zero()), 
    tloop(std::chrono::milliseconds::zero()),
    tsync(std::chrono::milliseconds::zero()),
    tpost_step_custom(std::chrono::milliseconds::zero()),
    tpost_step_base(std::chrono::milliseconds::zero())

  {
    int nx = this->mem->grid_size[0].length();
    int ny = this->mem->grid_size[1].length();
    int nz = this->mem->grid_size[2].length();
    w_LS.resize(nx,ny,nz);
    hgt_fctr_sclr.resize(nx,ny,nz);
    hgt_fctr_vctr.resize(nx,ny,nz);
    k_i.resize(nx, ny);
    r_l = 0.;

    blitz::thirdIndex k;
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
