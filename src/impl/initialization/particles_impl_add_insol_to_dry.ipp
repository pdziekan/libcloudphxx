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
    void particles_t<real_t, device>::impl::add_insol_to_dry(
      real_t rd_insol
    )
    {
      // add the insoluble aerosol
      const real_t rd3_insol = rd_insol * rd_insol * rd_insol;
      thrust::transform(
        rd3.begin()+n_part_old,
        rd3.end(),
        thrust::make_constant_iterator<real_t>(rd3_insol),
        thrust::plus<real_t>()
      );
    }
  };
};
