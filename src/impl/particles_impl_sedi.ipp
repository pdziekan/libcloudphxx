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
      template <class real_t>
      struct adve_euler_with_subs
      {
        const real_t dt, div_LS, z0;
        adve_euler_with_subs(const opts_init_t<real_t> &o) : dt(o.dt), div_LS(o.div_LS), z0(o.z0) {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &z, const thrust::tuple<const real_t, const real_t> tpl) // tpl: 0 - vt 1 - stretch
        {
          return z - dt / thrust::get<1>(tpl) * ( thrust::get<0>(tpl) + div_LS * (z - z0));  // Euler scheme (assuming vt positive!) NOTE!: we interpret here z0 as ground level, which might not be true! 
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sedi()
    {   
      namespace arg = thrust::placeholders;
 
      // settling due to sedimentation + large-scale subsidence = - div_LS * z
      thrust::transform(
        z.begin(), z.end(),                // input - 1st arg
        thrust::make_zip_iterator(   // input - 2nd arg
          thrust::make_tuple(
            vt.begin(),                    
            thrust::make_permutation_iterator(vert_stretch_prof.begin(), k.begin())
          )
        ),
        z.begin(),                         // output
        detail::adve_euler_with_subs<real_t>(opts_init)
      );
    }
  };  
};
