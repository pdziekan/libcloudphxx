// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct chem_vol_fun
      { // calculate drop volume
        const real_t pi;

        // ctor (pi() is not a __device__ function...)
        chem_vol_fun() :
          pi(boost::math::constants::pi<real_t>())
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw2) const
        {
#if !defined(__NVCC__)
	  using std::pow;
#endif
          return real_t(4./3) * pi * (pow(rw2, real_t(3./2)));
        }
      };

      template <class real_t>
      struct set_chem_flag
      { // set flag for those droplets that are big enough to have chemical reactions

        BOOST_GPU_ENABLED        
        unsigned int operator()(const real_t &rw2, const real_t &rd3)
        {
#if !defined(__NVCC__)
          using std::pow;
#endif
          return pow(rw2, real_t(3./2.)) > 1000. * rd3 ? 1 : 0;
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_vol_ante()
    {   
      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

      //calculate new drop volumes (to be used in Henrys law)
      thrust_device::vector<real_t> &V(tmp_device_real_part);

      thrust::transform(
        rw2.begin(), rw2.end(),         // input
        V.begin(),                      // output 
        detail::chem_vol_fun<real_t>()  // op
      );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_vol_post()
    {   
      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

      //save the current drop volume in V_old (to be used in the next step for Henrys law)
      thrust_device::vector<real_t> &V(tmp_device_real_part);
    
      thrust::copy(V.begin(), V.end(), V_old.begin());
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_flag_ante()
    { 
      thrust_device::vector<unsigned int> &chem_flag(tmp_device_n_part);
 
      // set flag to those SDs that take part in chemical reactions
      thrust::transform(
        rw2.begin(), rw2.end(),          // input 1st arg
        rd3.begin(),                     // input 2nd arg
        chem_flag.begin(),               // output
        detail::set_chem_flag<real_t>()  // op
      );
    }

  };  
};
