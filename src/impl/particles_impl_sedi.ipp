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
      template<class real_t>
      struct sedi_with_vt
      {
        const real_t dt;

        sedi_with_vt(real_t _dt): dt(_dt) {}

        BOOST_GPU_ENABLED
        real_t operator()(real_t z, real_t v) const
        {
          return z - dt * v;
        };
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sedi()
    {   
      namespace arg = thrust::placeholders;
 
      // settling due to sedimentation + large-scale subsidence
      thrust::transform(
        z.begin(), z.end(),                    // position
        vt.begin(),                                                    // terminal velocity 
        z.begin(),                         // output
        detail::sedi_with_vt<real_t>(opts_init.dt)
      );
    }
  };  
};
