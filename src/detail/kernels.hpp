#include <libcloudph++/common/moist_air.hpp>

#if defined(__NVCC__)
#  include <math_constants.h>
#endif

namespace libcloudphxx
{
  namespace lgrngn
  {

// global hall definitions...
     const std::array<float, 15> rad = {{6.e-6,8.e-6,10.e-6,15.e-6,20.e-6,25.e-6,30.e-6,40.e-6,50.e-6,60.e-6,70.e-6,100.e-6,150.e-6,200.e-6,300.e-6}};
     const std::array<float, 21> rat = {{0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5, 0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0}};

     const float eff[21][15] = 
     {
        {0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0},
        {0, 0, 0.0001, 0.0001, 0.0001,
        0, 0.000733, 0.001, 0.084, 0.05,
        0.2543, 0.5, 0.7807, 0.87, 0.97},
        {0, 0, 0.0001, 0.001733, 0.0001,
        0.00563, 0.002, 0.07, 0.4, 0.43,
        0.58, 0.79, 0.93, 0.96, 1},
        {0, 0, 0.014, 0.001733, 0.005,
        0.0156, 0.02667, 0.28, 0.62, 0.64,
        0.7629, 0.91, 0.9687, 0.98, 1},
        {0, 0, 0.014, 0.015, 0.016,
        0.028, 0.04, 0.5, 0.7, 0.77,
        0.84, 0.95, 0.95, 1, 1},
        {0, 0, 0.019, 0.02117, 0.022,
        0.0484, 0.1133, 0.62, 0.79, 0.84,
        0.8829, 0.95, 1, 1, 1},
        {0, 0, 0.019, 0.02983, 0.03,
        0.1226, 0.17, 0.68, 0.83, 0.87,
        0.9, 1, 1, 1, 1},
        {0, 0, 0.027, 0.02983, 0.043,
        0.1704, 0.3133, 0.74, 0.864, 0.89,
        0.9229, 1, 1, 1, 1},
        {0, 0, 0.027, 0.0395, 0.052,
        0.226, 0.4, 0.78, 0.88, 0.9,
        0.94, 1, 1, 1, 1},
        {0, 0, 0.033, 0.04883, 0.064,
        0.2708, 0.5167, 0.8, 0.9, 0.91,
        0.95, 1, 1, 1, 1},
        {0, 0, 0.033, 0.0555, 0.072,
        0.3184, 0.55, 0.8, 0.9, 0.91,
        0.95, 1, 1, 1, 1},
        {0, 0, 0.037, 0.0555, 0.079,
        0.3308, 0.5833, 0.8, 0.9, 0.91,
        0.95, 1, 1, 1, 1},
        {0, 0, 0.037, 0.0595, 0.082,
        0.336, 0.59, 0.78, 0.9, 0.91,
        0.95, 1, 1, 1, 1},
        {0, 0, 0.038, 0.05833, 0.08,
        0.3312, 0.5667, 0.77, 0.888, 0.91,
        0.95, 1, 1, 1, 1},
        {0, 0, 0.038, 0.05367, 0.076,
        0.3002, 0.54, 0.76, 0.88, 0.92,
        0.95, 1, 1, 1, 1},
        {0, 0, 0.036, 0.05367, 0.067,
        0.2855, 0.5033, 0.77, 0.882, 0.93,
        0.9743, 1, 1, 1, 1},
        {0, 0, 0.036, 0.0465, 0.057,
        0.2735, 0.49, 0.77, 0.89, 0.95,
        1, 1, 1, 1, 1},
        {0, 0, 0.032, 0.03967, 0.048,
        0.2619, 0.4633, 0.78, 0.938, 1,
        1.023, 1, 1, 1, 1},
        {0, 0, 0.032, 0.03267, 0.04,
        0.2476, 0.45, 0.79, 1.01, 1.03,
        1.04, 1, 1, 1, 1},
        {0, 0, 0.027, 0.03267, 0.033,
        0.2559, 0.4867, 0.95, 1.5, 1.7,
        2.543, 1, 1, 1, 1},
        {0, 0, 0.027, 0.027, 0.027,
        0.2735, 0.52, 1.4, 2.3, 3,
        4, 1, 1, 1, 1}
      };





    using detail::tpl_calc_wrap;

    template <typename real_t, typename n_t>
    struct kernel_base
    {
      // pointer to kernel parameters device vector
      thrust_device::pointer<real_t> k_params;

      // number of user-defined parameters
      n_t n_user_params; 

      // largest radius for which efficiency is defined, 0 - n/a
      real_t r_max;
   
      //ctor
      BOOST_GPU_ENABLED
      kernel_base(thrust_device::pointer<real_t> k_params, n_t n_user_params = 0, real_t r_max = 0.) : 
        k_params(k_params), n_user_params(n_user_params), r_max(r_max) {}

      // thrust requires that a default ctor exists
      kernel_base() = default;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &) const {return 0;}
    };


    //Golovin kernel
    template <typename real_t, typename n_t>
    struct kernel_golovin : kernel_base<real_t, n_t>
    {
      //ctor
      BOOST_GPU_ENABLED
      kernel_golovin(thrust_device::pointer<real_t> k_params) : kernel_base<real_t, n_t>(k_params, 1) {}

      // thrust requires that a default ctor exists
      kernel_golovin() = default;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };
#if !defined(__NVCC__)
        using std::max;
        using std::sqrt;
#endif

        const real_t rw2_a = thrust::get<rw2_a_ix>(tpl_wrap.get_rw());
        const real_t rw2_b = thrust::get<rw2_b_ix>(tpl_wrap.get_rw());

        real_t res =
#if !defined(__NVCC__)
        pi<real_t>()
#else
        CUDART_PI
#endif
        * 4. / 3.
        * kernel_base<real_t, n_t>::k_params[0]
        * max(
            thrust::get<n_a_ix>(tpl_wrap.get_rw()),
            thrust::get<n_b_ix>(tpl_wrap.get_rw())
          )
        * (
            rw2_a * sqrt(rw2_a) +
            rw2_b * sqrt(rw2_b) 
          );
        return res;
      }
    };


    //geometric kernel
    template <typename real_t, typename n_t>
    struct kernel_geometric : kernel_base<real_t, n_t>
    {
      //ctor (default one)
      BOOST_GPU_ENABLED
      kernel_geometric(thrust_device::pointer<real_t> k_params = thrust_device::pointer<real_t>(), n_t n_user_params = 0, real_t r_max = 0.) : 
        kernel_base<real_t, n_t>(k_params, n_user_params, r_max) {}

      //bilinear interpolation of collision efficiencies, required by dervied classes
      BOOST_GPU_ENABLED
      real_t interpolated_efficiency(real_t, real_t) const;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };
#if !defined(__NVCC__)
        using std::abs;
        using std::max;
        using std::sqrt;
#endif
        return 
#if !defined(__NVCC__)
        pi<real_t>()
#else
        CUDART_PI
#endif
        * max(
            thrust::get<n_a_ix>(tpl_wrap.get_rw()),
            thrust::get<n_b_ix>(tpl_wrap.get_rw())
          )
        * abs(
            thrust::get<vt_a_ix>(tpl_wrap.get_rw()) -
            thrust::get<vt_b_ix>(tpl_wrap.get_rw())
          )
        * (thrust::get<rw2_a_ix>(tpl_wrap.get_rw()) +
           thrust::get<rw2_b_ix>(tpl_wrap.get_rw()) +
           2.*sqrt(thrust::get<rw2_a_ix>(tpl_wrap.get_rw())*thrust::get<rw2_b_ix>(tpl_wrap.get_rw()))
          );
    //    return res;
      }
    };

    //geometric kernel with a multiplier
    template <typename real_t, typename n_t>
    struct kernel_geometric_with_multiplier : kernel_geometric<real_t, n_t>
    {
      //ctor
      BOOST_GPU_ENABLED
      kernel_geometric_with_multiplier(thrust_device::pointer<real_t> k_params) : kernel_geometric<real_t, n_t>(k_params, 1) {}

      // thrust requires that a default ctor exists
      kernel_geometric_with_multiplier() = default;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
        return kernel_geometric<real_t, n_t>::calc(tpl_wrap) * kernel_base<real_t, n_t>::k_params[0];
      }
    };

    //Long kernel
    template <typename real_t, typename n_t>
    struct kernel_long : kernel_geometric<real_t, n_t>
    {
      enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };
      //ctor
      BOOST_GPU_ENABLED
      kernel_long() : kernel_geometric<real_t, n_t>() {}

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
#if !defined(__NVCC__)
        using std::abs;
        using std::max;
        using std::min;
        using std::sqrt;
#endif
        real_t res = kernel_geometric<real_t, n_t>::calc(tpl_wrap);

        real_t r_L = max(sqrt(thrust::get<rw2_a_ix>(tpl_wrap.get_rw())), sqrt(thrust::get<rw2_b_ix>(tpl_wrap.get_rw())));
        if(r_L < 50.e-6)
        {
          real_t r_s = min(sqrt(thrust::get<rw2_a_ix>(tpl_wrap.get_rw())), sqrt(thrust::get<rw2_b_ix>(tpl_wrap.get_rw())));
          if(r_s <= 3e-6)
            res = 0.;
          else
            res *= 4.5e8 * r_L * r_L * (1. - 3e-6/r_s);
        }
        
        return  res;
      }
    };

    template <typename real_t, typename n_t>
    struct kernel_hall : kernel_geometric<real_t, n_t>
    {

      //ctor
      kernel_hall(thrust_device::pointer<real_t> k_params, real_t r_max) : kernel_geometric<real_t, n_t>(k_params, 0, r_max) {}

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };

#if !defined(__NVCC__)
        using std::sqrt;
#endif

        real_t radius = thrust::get<rw2_a_ix>(tpl_wrap.get_rw()) > thrust::get<rw2_b_ix>(tpl_wrap.get_rw()) ?
                     sqrt( thrust::get<rw2_a_ix>(tpl_wrap.get_rw()) ) : sqrt( thrust::get<rw2_b_ix>(tpl_wrap.get_rw()) ); 

        real_t ratio = thrust::get<rw2_a_ix>(tpl_wrap.get_rw()) > thrust::get<rw2_b_ix>(tpl_wrap.get_rw()) ?
                     sqrt( thrust::get<rw2_b_ix>(tpl_wrap.get_rw()) ) / sqrt( thrust::get<rw2_a_ix>(tpl_wrap.get_rw()) ) :
                     sqrt( thrust::get<rw2_a_ix>(tpl_wrap.get_rw()) ) / sqrt( thrust::get<rw2_b_ix>(tpl_wrap.get_rw()) );

        real_t efficiency; 
        // without interpolation
        /*
        int rad_pos, rat_pos;
        {
          auto lower_bound = std::lower_bound(rad.begin(), rad.end(), radius);
          if(lower_bound == rad.end()) rad_pos = 14;
          else if(lower_bound == rad.begin()) rad_pos = 0;
          else
          {
            if(*lower_bound - radius < radius - *(lower_bound-1)) rad_pos = lower_bound - rad.begin();
            else rad_pos = lower_bound - rad.begin() - 1;
          }
        }

        {
          auto lower_bound = std::lower_bound(rat.begin(), rat.end(), ratio);
          if(lower_bound == rat.end()) rat_pos = 20;
          else if(lower_bound == rat.begin()) rat_pos = 0;
          else
          {
            if(*lower_bound - ratio < ratio - *(lower_bound-1)) rat_pos = lower_bound - rat.begin();
            else rat_pos = lower_bound - rat.begin() - 1;
          }
        }
        efficiency = eff[rat_pos][rad_pos];
        */

        // with interpolation like in my_coad
        // find indexes of first not less value of radius and ratio
        int rad_pos = std::lower_bound(rad.begin(), rad.end(), radius) - rad.begin(); // [0, 15], could overflow
        int rat_pos = std::lower_bound(rat.begin(), rat.end(), ratio) - rat.begin(); // [0, 20]
        if(rad_pos < 15)
        {
          if(rad_pos >= 1)
          {
            real_t p = (radius - rad[rad_pos-1]) / (rad[rad_pos] - rad[rad_pos - 1]);
            real_t q = (ratio - rat[rat_pos-1]) / (rat[rat_pos] - rat[rat_pos - 1]);
            efficiency = (1-p)*(1-q)*eff[rat_pos-1][rad_pos-1] + 
                           p*(1-q)*eff[rat_pos-1][rad_pos] +
                           q*(1-p)*eff[rat_pos][rad_pos-1] + 
                           p*q*eff[rat_pos][rad_pos];
          }
          else
          {
            real_t q = (ratio - rat[rat_pos-1]) / (rat[rat_pos] - rat[rat_pos - 1]);
            efficiency = (1-q)*eff[rat_pos-1][0] + q*eff[rat_pos][0];
          }
        }
        else
        {
          real_t q = (ratio - rat[rat_pos-1]) / (rat[rat_pos] - rat[rat_pos - 1]);
          efficiency = (1-q)*eff[rat_pos-1][14] + q*eff[rat_pos][14];
          efficiency = efficiency > 1.? 1. : efficiency;
        }
/*

      do j=1,n r(j) - larger radius
      do i=1,j r(i) - smaller radius
         do k=2,15
            if (r(j).le.r0(k).and.r(j).ge.r0(k-1)) then
               ir=k
            elseif (r(j).gt.r0(15)) then
               ir=16
            elseif (r(j).lt.r0(1)) then
               ir=1
            endif
ir - indeks pierwszego wiekszego promienia od r(j)
         enddo
         rq=r(i)/r(j)
rq - stosunek malego do duzego
iq - indeks pierszego wiekszego stosunku
         do kk=2,21
            if (rq.le.rat(kk).and.rq.gt.rat(kk-1)) iq=kk
         enddo
         if (ir.lt.16) then         
            if (ir.ge.2) then
interpolacja wewnatrz
               p=(r(j)-r0(ir-1))/(r0(ir)-r0(ir-1))
               q=(rq-rat(iq-1))/(rat(iq)-rat(iq-1))
               ec(j,i)=(1.-p)*(1.-q)*ecoll(ir-1,iq-1)+
     &                 p*(1.-q)*ecoll(ir,iq-1)+
     &                 q*(1.-p)*ecoll(ir-1,iq)+
     &                 p*q*ecoll(ir,iq)
            else
ir==1, interpolacja tylko po rat
               q=(rq-rat(iq-1))/(rat(iq)-rat(iq-1))
               ec(j,i)=(1.-q)*ecoll(1,iq-1)+q*ecoll(1,iq)
            endif
         else
ir==16, interpolacja tylko po rat
            q=(rq-rat(iq-1))/(rat(iq)-rat(iq-1))
            ek=(1.-q)*ecoll(15,iq-1)+q*ecoll(15,iq)
            ec(j,i)=min(ek,1.d0)
         endif
         ec(i,j)=ec(j,i)
c         if (ec(i,j).lt.1.e-20) stop 99
      enddo
      enddo

*/
        return efficiency * kernel_geometric<real_t, n_t>::calc(tpl_wrap);
      }
    };

    template <typename real_t, typename n_t>
    struct kernel_geometric_with_efficiencies : kernel_geometric<real_t, n_t>
    {
      //ctor
      BOOST_GPU_ENABLED
      kernel_geometric_with_efficiencies(thrust_device::pointer<real_t> k_params, real_t r_max) : kernel_geometric<real_t, n_t>(k_params, 0, r_max) {}

      // thrust requires that a default ctor exists
      kernel_geometric_with_efficiencies() = default;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };

#if !defined(__NVCC__)
        using std::sqrt;
#endif

        return  kernel_geometric<real_t, n_t>::interpolated_efficiency(
                  sqrt( thrust::get<rw2_a_ix>(tpl_wrap.get_rw())),
                  sqrt( thrust::get<rw2_b_ix>(tpl_wrap.get_rw()))
                ) * kernel_geometric<real_t, n_t>::calc(tpl_wrap);
      }
    };

    // turbulent kernel from the 2015 JAS Onishi paper
    // two user params defined at initialization: 
    // turbulence dissipataion rate and Taylor microscale Reynolds number
    // TODO: get these values from flow characteristic (simulation-time during hskpng)
    //       cf. Benmoshe et al, JGR 2012
    template <typename real_t, typename n_t>
    struct kernel_onishi : kernel_geometric<real_t, n_t>
    {
      detail::wang_collision_enhancement_t<real_t> wang_collision_enhancement;

      //ctor
      BOOST_GPU_ENABLED
      kernel_onishi(thrust_device::pointer<real_t> k_params, real_t r_max) : kernel_geometric<real_t, n_t>(k_params, 1, r_max) {}

      // thrust requires that a default ctor exists
      kernel_onishi() = default;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };
        enum { rhod_ix, eta_ix, diss_rate_ix };

#if !defined(__NVCC__)
        using std::sqrt;
#endif
        real_t rwa = sqrt( thrust::get<rw2_a_ix>(tpl_wrap.get_rw()));
        real_t rwb = sqrt( thrust::get<rw2_b_ix>(tpl_wrap.get_rw()));
        real_t onishi_nograv = detail::kernel_onishi_nograv<real_t>           // value of the turbulent onishi kernel that neglects gravitational settling
        (
          rwa, rwb, kernel_base<real_t, n_t>::k_params[0], thrust::get<diss_rate_ix>(tpl_wrap.get_ro_calc()),                  // k_params[0] - Re_lambda
          thrust::get<eta_ix>(tpl_wrap.get_ro_calc()) / thrust::get<rhod_ix>(tpl_wrap.get_ro_calc()),                          // kinetic viscosity 
          common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres / thrust::get<rhod_ix>(tpl_wrap.get_ro_calc()) // ratio of water to air density
        );

        real_t geometric = kernel_geometric<real_t, n_t>::calc(tpl_wrap);
        real_t res = 
          kernel_geometric<real_t, n_t>::interpolated_efficiency(rwa, rwb) *             // stagnant air collision efficiency
          wang_collision_enhancement(rwa, rwb, kernel_base<real_t, n_t>::k_params[0]) *  // Wang turbulent collision efficiency enhancement, k_params[0] - epsilon
          sqrt(
            geometric * geometric +  // geometric kernel 
            onishi_nograv * onishi_nograv
          );

        return res;
      }
    };
  };
};


                            
