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
    p_0 = 101500 * si::pascals;
  const quantity<si::length, real_t> 
    z_0  = 0    * si::metres,
    Z    = 3000 * si::metres, // DYCOMS: 1500
    X    = 10000 * si::metres, // DYCOMS: 6400
    Y    = 10000 * si::metres; // DYCOMS: 6400
  const real_t z_i  = 1500; //initial inversion height
  const real_t heating_kappa = 85; // m^2/kg
  const real_t F_0 = 113; // w/m^2
  const real_t F_1 = 22; // w/m^2
  const real_t q_i = 8e-3; // kg/kg
  const real_t c_p = 1004; // J / kg / K
  const real_t z_abs = 2000; // [m] height above which absorber works

  const real_t D = 3.75e-6; // large-scale wind horizontal divergence [1/s]
  const real_t rho_i = 1.12; // kg/m^3

  const real_t F_sens = 15; //W/m^2, sensible heat flux
  const real_t F_lat = 115; //W/m^2, latent heat flux

  // liquid water potential temperature at height z
  template <class real_t>
  quantity<si::temperature, real_t> th_l(const real_t &z)
  {
    quantity<si::temperature, real_t> ret;
    ret = z < z_i ?
      289 * si::kelvins : 
      (303. + pow(z - z_i, 1./3)) * si::kelvins;
    return ret;
  }

  // water mixing ratio at height z
  struct r_t
  {
    quantity<si::dimensionless, real_t> operator()(const real_t &z) const
    {
      const quantity<si::dimensionless, real_t> q_t = z < z_i ?
        7.5e-3 : 0.5e-3; 
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

  // density profile as a function of altitude
  // hydrostatic and assuming constant theta (not used now)
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


  //aerosol lognormal dist. 
  const quantity<si::length, real_t>
    mean_rd = real_t(.08e-6) * si::metres;
  const quantity<si::dimensionless, real_t>
    sdev_rd = real_t(1.4);
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
    n_stp_clean = real_t(50e6) / si::cubic_metres,
    n_stp_polluted = real_t(300e6) / si::cubic_metres;

  //aerosol lognormal dist. for GCCN from Jorgen Jensen
  const quantity<si::length, real_t>
    mean_rd_gccn = real_t(.283e-6) * si::metres;
  const quantity<si::dimensionless, real_t>
    sdev_rd_gccn = real_t(2.235);
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
    n_stp_gccn = real_t(2.216e6) / si::cubic_metres;

  //aerosol chemical composition parameters (needed for activation)
  // for lgrngn:
  const quantity<si::dimensionless, real_t> kappa = .61; // CCN-derived value from Table 1 in Petters and Kreidenweis 2007
  // for blk_2m:
  const quantity<si::dimensionless, real_t> chem_b = .55; //ammonium sulphate //chem_b = 1.33; // sodium chloride

  //th, rv and surface fluxes relaxation time and height
  const quantity<si::time, real_t> tau_rlx = 1800 * si::seconds;
  const quantity<si::length, real_t> z_rlx_sclr = 200 * si::metres;


  // function expecting a libmpdata solver parameters struct as argument
  template <class T>
  void setopts(T &params, int nx, int nz)
  {
    params.dx = (X / si::metres) / (nx-1); 
    params.dz = (Z / si::metres) / (nz-1);
    params.di = params.dx;
    params.dj = params.dz;
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
      nz = solver.advectee().extent(z); 
    real_t 
      dx = (X / si::metres) / (nx-1), 
      dz = (Z / si::metres) / (nz-1); 

    // initial potential temperature & water vapour mixing ratio profiles
    solver.advectee(ix::th) = th_dry_fctr()(k * dz); 
    // randomly prtrb tht
    {
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
    }
    solver.advectee(ix::rv) = r_t()(k * dz); 
    // rndmly prtrb rv
    {
      std::random_device rd;
      auto seed = rd();
      std::mt19937 gen(seed);
      std::uniform_real_distribution<> dis(-0.2e-3, 0.2e-3);

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

      solver.advectee(ix::rv)(i_r, k_r) += prtrb(i_r, k_r);
    }

    solver.advectee(ix::u) = 0;
    solver.advectee(ix::w) = 0;  
   
    // absorbers
    solver.vab_coefficient() = where(k * dz >= z_abs,  1. / 100 * pow(sin(3.1419 / 2. * (k * dz - z_abs)/ (Z / si::metres - z_abs)), 2), 0);

    solver.vab_relaxed_state(0) = 0;
    solver.vab_relaxed_state(1) = 0;

    // density profile
    solver.g_factor() = rhod_fctr()(k * dz); 
  }

  // polluted lognormal aerosol distribution
  template <typename T>
  struct log_dry_radii_polluted : public libcloudphxx::common::unary_function<T>
  {
    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(mean_rd, sdev_rd, n_stp_polluted, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }

    log_dry_radii_polluted *do_clone() const 
    { return new log_dry_radii_polluted( *this ); }
  };

  // clean lognormal aerosol distribution
  template <typename T>
  struct log_dry_radii_clean : public libcloudphxx::common::unary_function<T>
  {
    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(mean_rd, sdev_rd, n_stp_clean, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }

    log_dry_radii_clean *do_clone() const 
    { return new log_dry_radii_clean( *this ); }
  };

  // polluted lognormal aerosol distribution with GCCN
  template <typename T>
  struct log_dry_radii_polluted_gccn : public libcloudphxx::common::unary_function<T>
  {
    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(mean_rd, sdev_rd, n_stp_polluted, quantity<si::dimensionless, real_t>(lnrd)) +
          lognormal::n_e(mean_rd_gccn, sdev_rd_gccn, n_stp_gccn, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }

    log_dry_radii_polluted_gccn *do_clone() const 
    { return new log_dry_radii_polluted_gccn( *this ); }
  };

  // clean lognormal aerosol distribution with GCCN
  template <typename T>
  struct log_dry_radii_clean_gccn : public libcloudphxx::common::unary_function<T>
  {
    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(mean_rd, sdev_rd, n_stp_clean, quantity<si::dimensionless, real_t>(lnrd)) +
          lognormal::n_e(mean_rd_gccn, sdev_rd_gccn, n_stp_gccn, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }

    log_dry_radii_clean_gccn *do_clone() const 
    { return new log_dry_radii_clean_gccn( *this ); }
  };
};
