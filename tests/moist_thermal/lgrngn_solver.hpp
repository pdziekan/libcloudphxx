#pragma once
#include <boost/assign/ptr_map_inserter.hpp> 

#include <libmpdata++/solvers/boussinesq.hpp>
#include <libmpdata++/output/gnuplot.hpp>
#include <libcloudph++/lgrngn/factory.hpp>
#include <future>
#include "setup.hpp"

using namespace libmpdataxx; // TODO: get rid of it?

template <class ct_params_t>
class lgrngn_solver : public 
  output::gnuplot<solvers::boussinesq<ct_params_t>>
{
  using parent_t = output::gnuplot<solvers::boussinesq<ct_params_t>>;

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

  // methods
  libcloudphxx::lgrngn::arrinfo_t<real_t> make_arrinfo(
    typename parent_t::arr_t arr
  ) {
    return libcloudphxx::lgrngn::arrinfo_t<real_t>(
      arr.dataZero(),
      arr.stride().data()
    );
  }


  void hook_ante_loop(int nt) 
  {
    parent_t::hook_ante_loop(nt); 
    // init particles
    if (this->rank == 0)
    {
      assert(params.dt != 0);

      params.backend = libcloudphxx::lgrngn::OpenMP;

      params.cloudph_opts_init.dt = params.dt; // advection timestep = microphysics timestep
      params.cloudph_opts_init.dx = params.di;
      params.cloudph_opts_init.dz = params.dj;

      // libmpdata++'s grid interpretation
      params.cloudph_opts_init.x0 = params.di / 2;
      params.cloudph_opts_init.z0 = params.dj / 2;
      params.cloudph_opts_init.x1 = (this->mem->grid_size[0].length() - .5) * params.di;
      params.cloudph_opts_init.z1 = (this->mem->grid_size[1].length() - .5) * params.dj;
      params.cloudph_opts_init.nx = (this->mem->grid_size[0].length());
      params.cloudph_opts_init.nz = (this->mem->grid_size[1].length());

      params.cloudph_opts_init.sd_conc = 64;
      params.cloudph_opts_init.n_sd_max = 
        params.cloudph_opts_init.sd_conc * 
        params.cloudph_opts_init.nx * 
        params.cloudph_opts_init.nz; 

      params.cloudph_opts_init.terminal_velocity = libcloudphxx::lgrngn::vt_t::khvorostyanov_spherical;
      params.cloudph_opts_init.kernel = libcloudphxx::lgrngn::kernel_t::hall_davis_no_waals;

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

      prtcls->init(
        make_arrinfo(this->mem->advectee(ix::tht)),
        make_arrinfo(this->mem->advectee(ix::rv)),
       // make_arrinfo(this->mem->g_factor())
      );
    }
  }

  void hook_ante_step()
  {
    parent_t::hook_ante_step(); 
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
  }  
};
