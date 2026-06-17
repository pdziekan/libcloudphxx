// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_kappa(
      const real_t kappa,     // kappa of the soluble part
      const real_t rd_insol  // size of the insoluble part
    )
    {
      // assuming that rd3 = rd3_sol + rd3_insol is already filled

      const real_t rd3_insol = rd_insol * rd_insol * rd_insol;
      
      thrust::transform(
        rd3.begin() + n_part_old, rd3.end(),
        kappa.begin(),
        (arg::_1 - rd3_insol) / arg::_1 * kappa
        // assuming that soluble and insoluble densities are the same (hence using volume fraction mixing, not mass fraction as should (?) be done)
        
    }
  };
};

