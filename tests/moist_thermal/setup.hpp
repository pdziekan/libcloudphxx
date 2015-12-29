#pragma once

#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>

namespace setup
{
  using real_t = float;

  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;
  namespace moist_air = libcloudphxx::common::moist_air;
  namespace const_cp = libcloudphxx::common::const_cp;

  quantity<si::dimensionless, real_t> RH_to_rv(const real_t &RH, const quantity<si::temperature, real_t> &T, const quantity<si::pressure, real_t> &p)
  {
    return moist_air::eps<real_t>() * RH * const_cp::p_vs<real_t>(T) / (p - RH * const_cp::p_vs<real_t>(T));
  }

  real_t env_RH = 0.2;
  real_t prtrb_RH = 1.0200;

  const quantity<si::temperature, real_t>
    T_0(283. * si::kelvins);  // surface temperature
  const quantity<si::pressure, real_t>
    p_0 = 85000 * si::pascals; // surface pressure
  const quantity<si::dimensionless, real_t> rv_0(RH_to_rv(env_RH, T_0, p_0));
  // theta (std) at surface
  const quantity<si::temperature, real_t> th_0 = T_0 * pow(100000 * si::pascals / p_0,  moist_air::R_d<real_t>() / moist_air::c_pd<real_t>() );
  const quantity<si::temperature, real_t> th_0_dry = theta_dry::std2dry<real_t>(th_0, rv_0);
  const quantity<si::length, real_t> z_0(0. * si::metres);

  // calc theta std at height z (c.f. Grabowski Clark 1991)
  template <class real_t>
  quantity<si::temperature, real_t> th(const real_t &z)
  {
    // theta at height z
    const quantity<si::temperature, real_t> th = ((th_0 / si::kelvins * real_t(exp(1.3e-5 * z))) * si::kelvins);
    return th;
  }

  // assume hydrostatic pressure
  template <class real_t>
  quantity<si::pressure, real_t> p(const real_t &z)
  {
    quantity<si::pressure, real_t> p = hydrostatic::p(
      z * si::metres, th_0_dry, rv_0, z_0, p_0
    ); // aprox; it should be anelastic; also theta and rv are not constant
    return p;
  }

  // theta std to temperature
  template <class real_t>
  quantity<si::temperature, real_t> th2T(const quantity<si::temperature, real_t> &th, const quantity<si::pressure, real_t> &p)
  {
    quantity<si::temperature, real_t> T = th  * pow(100000 * si::pascals / p,  - moist_air::R_d<real_t>() / moist_air::c_pd<real_t>() );
    return T;
  }

  struct env_rv
  {
    quantity<si::dimensionless, real_t> operator()(const real_t &z) const
    {
      return RH_to_rv(env_RH, th2T(th(z), p(z)), p(z));
    }
  BZ_DECLARE_FUNCTOR(env_rv);
  };

  struct prtrb_rv
  {
    quantity<si::dimensionless, real_t> operator()(const real_t &z) const
    {
      return RH_to_rv(prtrb_RH, th2T(th(z), p(z)), p(z));
    }
  BZ_DECLARE_FUNCTOR(prtrb_rv);
  };

  struct th_dry
  {
    real_t operator()(const real_t &z) const
    {
      return theta_dry::std2dry<real_t>(th(z), env_rv()(z)) / si::kelvins;
    }
  BZ_DECLARE_FUNCTOR(th_dry);
  };

  // density profile as a function of altitude
  // also approx.
  struct rhod
  {
    real_t operator()(real_t z) const
    {
      quantity<si::mass_density, real_t> rhod = theta_std::rhod(
        p(z), th(z), env_rv()(z)
      );

      return rhod / si::kilograms * si::cubic_metres;
    }

    // to make the rhod() functor accept Blitz arrays as arguments
    BZ_DECLARE_FUNCTOR(rhod);
  };

  //aerosol bimodal lognormal dist. 
  const quantity<si::length, real_t>
    mean_rd1 = real_t(.04e-6 / 2) * si::metres,
    mean_rd2 = real_t(.15e-6 / 2) * si::metres;
  const quantity<si::dimensionless, real_t>
    sdev_rd1 = real_t(1.4),
    sdev_rd2 = real_t(1.6);
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
    n1_stp = real_t(60e6) / si::cubic_metres,
    n2_stp = real_t(40e6) / si::cubic_metres;
  const quantity<si::dimensionless, real_t> kappa = .61; // CCN-derived value from Table 1 in Petters and Kreidenweis 2007


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

