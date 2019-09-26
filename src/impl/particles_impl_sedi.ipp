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
        arg::_1 - opts_init.dt * (thrust::get<0>(arg::_2) / thrust::get<1>(arg::_2) + opts_init.div_LS * (arg::_1 - opts_init.z0))   // Euler scheme (assuming vt positive!) NOTE!: we interpret here z0 as ground level, which might not be true! 
      );
    }
  };  
};
