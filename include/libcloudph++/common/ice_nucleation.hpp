#pragma once

#include "units.hpp"
#include "const_cp.hpp"
#include <algorithm>
#include <thrust/tuple.h>

#if defined(__NVCC__)
#  include <math_constants.h>
#endif

namespace libcloudphxx
{
  namespace common
  {
    namespace ice_nucleation
    {
      enum class INP_t {mineral, AgI}; // types of ice nucleating particles, TODO: add more types

      // Inverse CDF for singular freezing temperature as defined in eq. 1 in Shima et al., 2020
      // frozen_fraction(T) = p(T_f > T) 
      // CDF of T_f = 1 - frozen_fraction(T)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::temperature, real_t> T_freeze_CDF_inv(
      const INP_t INP_type,      // type of ice nucleating particle
      const real_t rd2_insol,     // radius squared of insoluble particle in m^2
      const real_t rand           // random number between [0, 1]
        ) {
        static constexpr quantity<si::temperature, real_t> T_homo = real_t(235.15) * si::kelvin;
        // static constexpr quantity<si::temperature, real_t> Niemand_T_min = real_t(273.15 - 36) * si::kelvin;
        // static constexpr quantity<si::temperature, real_t> Niemand_T_max = real_t(273.15 - 12) * si::kelvin;
        static constexpr quantity<si::dimensionless, real_t> Omanovic_b = 0.97;
        static constexpr quantity<si::dimensionless, real_t> Omanovic_k = 0.88;
        static constexpr quantity<si::temperature, real_t> Omanovic_T0 = real_t(263.95) * si::kelvin;

        if(rd2_insol < 1e-20)  return T_homo;

        switch(INP_type)
        {
          case INP_t::mineral: {
            // active site ns(T) parameterization from Niemand et al. 2012 for mineral dust
            // Shima et al. 2020 (and many others): p(T_f > T) = 1 - exp(A * ns(T)) 
            // Niemand et al. 2012 for mineral dust: ns(T) = exp(-0.517(T - 273.15) + 8.934) [m^-2]
            // NOTE (not used by us): Shima et al. 2020: use Niemand only in the range -36 C to -12 C, above ns(T)=0 and below ns(T)=ns(-36)

            const real_t A = real_t(4)
            #if !defined(__NVCC__)
                * pi<real_t>()
            #else 
                * CUDART_PI
            #endif
            * rd2_insol; // surface area of the insoluble particle

            const real_t Niemand_T_freeze = real_t(real_t(273.15) + (real_t(8.934) - log(- log(real_t(1.) - rand) / (real_t(4) * pi<real_t>() * rd2_insol)) ) / real_t(0.517));
            return std::max(real_t(T_homo / si::kelvin), Niemand_T_freeze) * si::kelvin;
            break;
          }
          case INP_t::AgI: {
            // Omanovic et al. 2024 for AgI: frozen_fraction = b [ 1 - 1 / (1 + exp(-k(T-T0)))]; b=0.97, k=0.88, T0=263.95K ; applicable for rd_insol>20nm
            // leads to T_f = T0 - 1/k ln((1-R)/(b-1+R)) , where R is uniformly distributed [0,1]; as T->0 FF=b, so if R>b, set Tf=0?
            if(rand > Omanovic_b) return T_homo;
            const real_t Omanovic_T_freeze = real_t(Omanovic_T0/ si::kelvin) - (real_t(1.) / Omanovic_k) * log((real_t(1.) - rand) / (Omanovic_b - real_t(1.) + rand));
            return std::max(real_t(T_homo / si::kelvin), Omanovic_T_freeze) * si::kelvin;
            break;
          }
          default:
            throw std::runtime_error("Unrecognized INP type");
        }
      }


      template<typename real_t>
      struct T_freeze_CDF_inv_functor
      {
        INP_t INP_type;

        T_freeze_CDF_inv_functor(INP_t INP_type)
          : INP_type(INP_type) {}

        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t> &tpl) const
        {
          const real_t &rd2_insol = thrust::get<0>(tpl);  // from rd2 vector
          const real_t &rand         = thrust::get<1>(tpl);  // from rand vector

          return ice_nucleation::template T_freeze_CDF_inv<real_t>(
            INP_type,
            rd2_insol,
            rand
          ) / si::kelvin;
        }
      };

      // Probability of time-dependent freezing,
      // heterogeneous as in Arabas et al., 2025 and homogeneous as in Koop & Murray, 2016
      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t p_freeze(
      const INP_t& INP_type,     // type of ice nucleating particle
      const real_t rd2_insol,    // radius squared of insoluble particle in m^2
      const real_t rw2,          // wet radius squared in m^2
      const real_t T,            // temperature in kelvin
      const real_t dt            // time step in seconds
        )
      {
        if (rd2_insol > real_t(0))
        {
          real_t A = real_t(4)
          #if !defined(__NVCC__)
              * pi<real_t>()
          #else
              * CUDART_PI
          #endif
          * rd2_insol; // surface area of the insoluble particle
          real_t d_aw = real_t(1) - const_cp::p_vsi<real_t>(T * si::kelvin)/ const_cp::p_vs<real_t>(T * si::kelvin); // water activity
          if (INP_type == INP_t::mineral)
          {
            real_t J_het = pow(real_t(10), real_t(-1.35) + real_t(22.62) * d_aw) * real_t(1e4); // nucleation rate
            return 1 - exp(- J_het * A * dt);
          }
          else if (INP_type == INP_t::AgI)
          {
            throw std::runtime_error("AgI time-dependent freezing not implemented yet");
            return 0;
          }
          else
          {
            throw std::runtime_error("Unrecognized INP type");
            return 0;
          }
        }
        else
        {
          real_t V = real_t(4./3.)
          #if !defined(__NVCC__)
              * pi<real_t>()
          #else
              * CUDART_PI
          #endif
          * pow(rw2, real_t(3./2.)); // droplet volume

          real_t dT = T - real_t(273.15);
          real_t x = - real_t(3020.684) - real_t(425.921)*pow(dT,real_t(1)) - real_t(25.9779)*pow(dT,real_t(2))
                      - real_t(0.868451)*pow(dT,real_t(3)) - real_t(0.0166203)*pow(dT,real_t(4))
                      - real_t(0.000171736)*pow(dT,real_t(5)) - real_t(0.000000746953)*pow(dT,real_t(6));
          real_t J_hom = pow(real_t(10), x) * real_t(1e6); // nucleation rate
          return 1 - exp(- J_hom * V * dt);
        }
      }


      template<typename real_t>
      struct p_freeze_functor
      {
        INP_t INP_type;
        real_t dt;

        p_freeze_functor(INP_t INP_type, real_t dt)
          : INP_type(INP_type), dt(dt) {}

        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) const
        {
          const real_t &rd2_insol = thrust::get<0>(tpl);  // radius squared of insoluble particle
          const real_t &rw2       = thrust::get<1>(tpl);  // wet radius squared
          const real_t &T         = thrust::get<2>(tpl);  // temperature in kelvin

          return ice_nucleation::p_freeze<real_t>(
            INP_type,
            rd2_insol,
            rw2,
            T,
            dt
          );
        }
      };

    };
  };
};
