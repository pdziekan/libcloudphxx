// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <thrust/sequence.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    // init number of SDs to be initialized per cell
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_count_num()
    {
      if(n_dims > 0)
      {
        namespace arg = thrust::placeholders;
        // some cells may be used only partially in the Lagrangian method;
        // sd_conc defines number of SDs per Eulerian cell
        thrust::transform(dv.begin(), dv.end(), count_num.begin(), (real_t(opts_init.sd_conc) * arg::_1 / (opts_init.dx * opts_init.dy * opts_init.dz) + real_t(0.5))); 
      }
      // parcel setup
      else
        thrust::fill(count_num.begin(), count_num.end(), opts_init.sd_conc);
    }
  };
};
