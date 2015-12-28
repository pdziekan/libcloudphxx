#pragma once
#include <boost/assign/ptr_map_inserter.hpp> 

#include <libmpdata++/solvers/mpdata_rhs_vip_prs.hpp>
#include <libmpdata++/output/hdf5.hpp>
#include <libcloudph++/lgrngn/factory.hpp>
#include <future>
#include "setup.hpp"

using namespace libmpdataxx; // TODO: get rid of it?

template <class ct_params_t>
class lgrngn_solver : public 
  output::hdf5<solvers::mpdata_rhs_vip_prs<ct_params_t>>
{
  using parent_t = output::hdf5<solvers::mpdata_rhs_vip_prs<ct_params_t>>;

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;

  struct rt_params_t : parent_t::rt_params_t 
  { 
    int backend;
    libcloudphxx::lgrngn::opts_t<real_t> cloudph_opts;
    libcloudphxx::lgrngn::opts_init_t<real_t> cloudph_opts_init;
  };

  // member fields
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto_t<real_t>> prtcls;

  rt_params_t params;
  blitz::Array<real_t,2> rhod;
  blitz::Array<real_t,2> tht_env_init;
//  std::future<real_t> ftr;

  // methods
  libcloudphxx::lgrngn::arrinfo_t<real_t> make_arrinfo(
    typename parent_t::arr_t arr
  ) {
    return libcloudphxx::lgrngn::arrinfo_t<real_t>(
      arr.dataZero(),
      arr.stride().data()
    );
  }

  // helper methods
  void diag()
  {
    assert(this->rank == 0);

    // recording super-droplet concentration per grid cell 
    prtcls->diag_sd_conc();
    this->record_aux("sd_conc", prtcls->outbuf());

    // cloud water
    prtcls->diag_wet_rng(.5e-6, 25e-6);
    prtcls->diag_wet_mom(3);
    this->record_aux("rl_mom3", prtcls->outbuf());

    // rain water
    prtcls->diag_wet_rng(25e-6, 1);
    prtcls->diag_wet_mom(3);
    this->record_aux("rr_mom3", prtcls->outbuf());
  } 

  void update_rhs(
    libmpdataxx::arrvec_t<
      typename parent_t::arr_t
    > &rhs, 
    const real_t &dt, 
    const int &at 
  ) {
    parent_t::update_rhs(rhs, dt, at); 

    const auto &Tht = this->state(ix::tht); 
    const auto &ijk = this->ijk;

    if(at>1)
    {
      std::cout << "update at" << at << " tht init" << tht_env_init << std::flush;
      rhs.at(ix::w)(ijk) += 
        9.81 * (Tht(ijk) - tht_env_init(ijk)) / tht_env_init(ijk); 
    }
  }

  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 
    // init particles
    if (this->rank == 0)
    {
      assert(params.dt != 0);

      params.backend = libcloudphxx::lgrngn::OpenMP;
      params.cloudph_opts.cond = true;
      params.cloudph_opts.adve = true;
      params.cloudph_opts.coal = true;
      params.cloudph_opts.sedi = true;

      params.cloudph_opts_init.dt = params.dt; // advection timestep = microphysics timestep
      params.cloudph_opts_init.dx = params.di;

      // libmpdata++'s grid interpretation
      params.cloudph_opts_init.x0 = params.di / 2;
      params.cloudph_opts_init.z0 = params.dj / 2;
      params.cloudph_opts_init.x1 = (this->mem->grid_size[0].length() - .5) * params.di;
      params.cloudph_opts_init.z1 = (this->mem->grid_size[1].length() - .5) * params.dj;

      params.cloudph_opts_init.sstp_coal = 1;
      params.cloudph_opts_init.sstp_cond = 10;

      params.cloudph_opts_init.sd_conc = 64;
      params.cloudph_opts_init.n_sd_max = 
        params.cloudph_opts_init.sd_conc * 
        params.cloudph_opts_init.nx * 
        params.cloudph_opts_init.nz; 

      params.cloudph_opts_init.terminal_velocity = libcloudphxx::lgrngn::vt_t::khvorostyanov_nonspherical;
      params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::geometric;

      boost::assign::ptr_map_insert<
        setup::log_dry_radii<setup::real_t> // value type
      >(
        params.cloudph_opts_init.dry_distros // map
      )(
        setup::kappa // key
      );

      prtcls.reset(libcloudphxx::lgrngn::factory<real_t>(
        (libcloudphxx::lgrngn::backend_t)params.backend,
        params.cloudph_opts_init
      ));

      rhod.resize(params.cloudph_opts_init.nx, params.cloudph_opts_init.nz); 

      blitz::secondIndex j;
      rhod = setup::rhod()(j * params.cloudph_opts_init.dz);

      prtcls->init(
        make_arrinfo(this->mem->advectee(ix::tht)),
        make_arrinfo(this->mem->advectee(ix::rv)),
        make_arrinfo(rhod)
      );
      diag();
    }
  }

  void hook_ante_step()
  {
    parent_t::hook_ante_step(); 

    this->mem->barrier();

  }

  void hook_post_step()
  {
    parent_t::hook_post_step(); // includes output

    this->mem->barrier();

    if (this->rank == 0)
    {
      // assuring previous async step finished ...
/*      if (
        this->timestep != 0 && // ... but not in first timestep ...
        ((this->timestep - 1) % this->outfreq != 0) // ... and not after diag call
      ) {
        assert(ftr.valid());
        ftr.get();
      } 
*/
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
       /* {
          blitz::secondIndex j;
          Cx /= rhod()(   j     * this->dj);
          Cz /= rhod()((j - .5) * this->dj);
        }*/

//        std::cout << this->mem->advectee(ix::rv);
        // running synchronous stuff
        prtcls->step_sync(
          params.cloudph_opts,
          make_arrinfo(this->mem->advectee(ix::tht)),
          make_arrinfo(this->mem->advectee(ix::rv)),
          make_arrinfo(rhod),
//          make_arrinfo(this->mem->advectee(ix::u)),
          make_arrinfo(Cx), // ix::u ?
          libcloudphxx::lgrngn::arrinfo_t<real_t>(),
  //        make_arrinfo(this->mem->advectee(ix::w))
          make_arrinfo(Cz) // ix:w ?
        );
//        std::cout << this->mem->advectee(ix::rv);
      }

      // running asynchronous stuff
      { 
        using libcloudphxx::lgrngn::particles_t;
        using libcloudphxx::lgrngn::OpenMP;

        prtcls->step_async(params.cloudph_opts);
//        assert(!ftr.valid());
    /*    ftr = std::async(
          std::launch::async,  
          &particles_t<real_t, OpenMP>::step_async, 
          dynamic_cast<particles_t<real_t, OpenMP>*>(prtcls.get()),
          params.cloudph_opts
        );*/
//        assert(ftr.valid());
      }
      
      // performing diagnostics
      if (this->timestep % this->outfreq == 0)
      { 
//        assert(ftr.valid());
  //      ftr.get();
        diag();
      }
    }
  }

  public:

  // ctor
  lgrngn_solver( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    params(p)
  {
    params.cloudph_opts_init.nx = (this->mem->grid_size[0].length());
    params.cloudph_opts_init.nz = (this->mem->grid_size[1].length());
    params.cloudph_opts_init.dz = params.dj;
    tht_env_init.resize(params.cloudph_opts_init.nx, params.cloudph_opts_init.nz); 
    blitz::secondIndex j;
    tht_env_init = setup::th_dry()(j * params.cloudph_opts_init.dz);

    std::cout << "ctr tht init" << tht_env_init << std::flush;
  }  
};
