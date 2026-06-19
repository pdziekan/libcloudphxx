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
    void particles_t<real_t, device>::impl::init_insol(
      real_t soluble_fraction // volume fraction of the soluble part
    )
    {
      // rd2_insol[i] = ((1 - soluble_fraction) * rd3[i])^(2/3)
      namespace arg = thrust::placeholders;
      thrust::transform(
        rd3.begin() + n_part_old,
        rd3.end(),
        rd2_insol.begin() + n_part_old,
        [soluble_fraction] BOOST_GPU_ENABLED (real_t rd3_val) {
          return std::pow((real_t(1) - soluble_fraction) * rd3_val, real_t(2) / real_t(3));
        }
      );
    }
  };
};
