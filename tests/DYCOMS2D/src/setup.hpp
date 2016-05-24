#pragma once

#include <iostream>

#include <blitz/array.h> 

#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>

// TODO: relaxation terms still missing

// 8th ICMW case 1 by Wojciech Grabowski)
namespace setup 
{
  using real_t = float;

  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;

  enum {x, z}; // dimensions
  const quantity<si::pressure, real_t> 
    p_0 = 101780 * si::pascals;
  const quantity<si::length, real_t> 
    z_0  = 0    * si::metres,
    Z    = 1500 * si::metres, // DYCOMS: 1500
    X    = 6400 * si::metres, // DYCOMS: 6400
    Y    = 6400 * si::metres; // DYCOMS: 6400
  const real_t z_i  = 795; //initial inversion height
  const real_t heating_kappa = 85; // m^2/kg
  const real_t F_0 = 70; // w/m^2
  const real_t F_1 = 22; // w/m^2
  const real_t q_i = 8e-3; // kg/kg
  const real_t c_p = 1004; // J / kg / K
  const real_t z_abs = 1200; // [m] height above which absorber works

  const real_t D = 3.75e-6; // large-scale wind horizontal divergence [1/s]
  const real_t rho_i = 1.12; // kg/m^3

  const real_t F_sens = 16; //W/m^2, sensible heat flux
  const real_t F_lat = 93; //W/m^2, latent heat flux
  const real_t u_fric = 0.25; // m/s, friction velocity

  // liquid water potential temperature at height z
  template <class real_t>
  quantity<si::temperature, real_t> th_l(const real_t &z)
  {
    quantity<si::temperature, real_t> ret;
    ret = z < z_i ?
      288.3 * si::kelvins : 
      (295. + pow(z - z_i, 1./3)) * si::kelvins;
    return ret;
  }

  // water mixing ratio at height z
  struct r_t
  {
    quantity<si::dimensionless, real_t> operator()(const real_t &z) const
    {
      const quantity<si::dimensionless, real_t> q_t = z < z_i ?
        9.45e-3 : 
        (5. - 3. * (1. - exp((z_i - z)/500.))) * 1e-3;
      return q_t;
    }
    BZ_DECLARE_FUNCTOR(r_t);
  };

  // initial dry air potential temp at height z, assuming theta_std = theta_l (spinup needed)
  struct th_dry_fctr
  {
    real_t operator()(const real_t &z) const
    {
      return theta_dry::std2dry<real_t>(th_l(z), r_t()(z)) / si::kelvins;
    }
    BZ_DECLARE_FUNCTOR(th_dry_fctr);
  };

  // westerly wind
  struct u
  {
    real_t operator()(const real_t &z) const
    {
      return 3. + 4.3 * z / 1000.; 
    }
    BZ_DECLARE_FUNCTOR(u);
  };

  // southerly wind
  struct v
  {
    real_t operator()(const real_t &z) const
    {
      return -9. + 5.6 * z / 1000.; 
    }
    BZ_DECLARE_FUNCTOR(v);
  };

  // large-scale vertical wind
  struct w_LS_fctr
  {
    real_t operator()(const real_t &z) const
    {
      return -D * z; 
    }
    BZ_DECLARE_FUNCTOR(w_LS_fctr);
  };

  // density profile as a function of altitude
  // hydrostatic
  struct rhod_fctr
  {
    real_t operator()(real_t z) const
    {
      quantity<si::pressure, real_t> p = hydrostatic::p(
	z * si::metres, th_dry_fctr()(0.) * si::kelvins, r_t()(0.), z_0, p_0
      );
      
      quantity<si::mass_density, real_t> rhod = theta_std::rhod(
	p, th_dry_fctr()(0.) * si::kelvins, r_t()(0.)
      );

      return rhod / si::kilograms * si::cubic_metres;
    }

    // to make the rhod() functor accept Blitz arrays as arguments
    BZ_DECLARE_FUNCTOR(rhod_fctr);
  };


  //aerosol bimodal lognormal dist. 
  const quantity<si::length, real_t>
    mean_rd1 = real_t(.011e-6) * si::metres,
    mean_rd2 = real_t(.06e-6) * si::metres;
  const quantity<si::dimensionless, real_t>
    sdev_rd1 = real_t(1.2),
    sdev_rd2 = real_t(1.7);
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
    n1_stp = real_t(125e6) / si::cubic_metres,
    n2_stp = real_t(65e6) / si::cubic_metres;

  //aerosol chemical composition parameters (needed for activation)
  // for lgrngn:
  const quantity<si::dimensionless, real_t> kappa = .61; // CCN-derived value from Table 1 in Petters and Kreidenweis 2007
  // for blk_2m:
  const quantity<si::dimensionless, real_t> chem_b = .55; //ammonium sulphate //chem_b = 1.33; // sodium chloride

  //th, rv and surface fluxes relaxation time and height
  const quantity<si::time, real_t> tau_rlx = 300 * si::seconds;
  const quantity<si::length, real_t> z_rlx = 50 * si::metres;

  // function expecting a libmpdata solver parameters struct as argument
  template <class T>
  void setopts(T &params, int nx, int ny, int nz)
  {
    params.dx = (X / si::metres) / (nx-1); 
    params.dz = (Z / si::metres) / (nz-1);
    params.di = params.dx;
    params.dj = params.dz;
 
    // prescribed density
/*    blitz::Array<real_t, 2> rhod(nx, nz);
    {
      blitz::secondIndex j;
      rhod = rhod_fctr()(j * params.dz);
    }
    std::cout << "rhod w setopts" << rhod << std::endl;

    params.rhod = new blitz::Array<real_t, 2>(rhod.dataFirst(), rhod.shape(), blitz::neverDeleteData);
    std::cout << "params.rhod w setopts" << params.rhod << " " << *params.rhod << std::endl;

    // prescribed large-scale vertical wind
    blitz::Array<real_t, 2> w_LS(nx, nz);
    {
      blitz::secondIndex j;
      w_LS = w_LS_fctr()(j * params.dz);
    }
    params.w_LS = new blitz::Array<real_t, 2>(w_LS.dataFirst(), w_LS.shape(), blitz::neverDeleteData);

    // prescribed initial temp profile
    blitz::Array<real_t, 2> th_init(nx, nz);
    {
      blitz::secondIndex j;
      th_init = th_dry_fctr()(j * params.dz);
    }
    params.th_init = new blitz::Array<real_t, 2>(th_init.dataFirst(), th_init.shape(), blitz::neverDeleteData);
*/
  }

  // function expecting a libmpdata++ solver as argument
  template <class concurr_t>
  void intcond(concurr_t &solver)
  {
    using ix = typename concurr_t::solver_t::ix;

    // helper ondex placeholders
    blitz::firstIndex i;
    blitz::secondIndex k;

    // dx, dy ensuring 1500x1500 domain
    int 
      nx = solver.advectee().extent(x), 
      nz = solver.advectee().extent(y); 
    real_t 
      dx = (X / si::metres) / (nx-1), 
      dz = (Z / si::metres) / (nz-1); 

    // initial potential temperature & water vapour mixing ratio profiles
    solver.advectee(ix::th) = th_dry_fctr()(k * dz); 
    // randomly prtrb tht
    std::random_device rd;
    auto seed = rd();
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(-0.1, 0.1);

    blitz::Array<real_t, 2> prtrb(nx, nz);
    for (int ii = 0; ii < nx; ++ii)
    {
      for (int kk = 0; kk < nz; ++kk)
      {
         prtrb(ii, kk) = dis(gen);
      }
    }
    auto i_r = blitz::Range(0, nx - 1);
    auto k_r = blitz::Range(0, nz - 1);

    // enforce cyclic perturbation
    prtrb(nx - 1, k_r) = prtrb(0, k_r);

    solver.advectee(ix::th)(i_r, k_r) += prtrb(i_r, k_r);

    solver.advectee(ix::rv) = r_t()(k * dz); 

    solver.advectee(ix::u) = 0;
    solver.advectee(ix::u)(i_r, k_r)= setup::u()(k * dz);
    solver.advectee(ix::w) = 0;  
   
    // absorbers
    solver.vab_coefficient() = where(k * dz >= z_abs, 1. / 1020 * (k * dz - z_abs) / (Z / si::metres - z_abs), 0);
    solver.vab_relaxed_state(0) = solver.advectee(ix::u);
    solver.vab_relaxed_state(1) = 0;

    // density profile
    solver.g_factor() = rhod_fctr()(k * dz); // TODO: reenable g_factor (and nug option) once it works in 3D libmpdata
  }

  // lognormal aerosol distribution
  template <typename T>
  struct log_dry_radii : public libcloudphxx::common::unary_function<T>
  {
    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(mean_rd1, sdev_rd1, n1_stp, quantity<si::dimensionless, real_t>(lnrd)) +
          lognormal::n_e(mean_rd2, sdev_rd2, n2_stp, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }

    log_dry_radii *do_clone() const 
    { return new log_dry_radii( *this ); }
  };
};
