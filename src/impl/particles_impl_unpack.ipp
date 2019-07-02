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
    void particles_t<real_t, device>::impl::unpack_n(const int &n_copied)
    {
      n_part_old = n_part;
      n_part += n_copied;
      assert(opts_init.n_sd_max >= n_part);
      n.resize(n_part);
      thrust::copy(in_n_bfr.begin(), in_n_bfr.begin() + n_copied, n.begin() + n_part_old);
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::unpack_real(const int &n_copied)
    {
      auto it = distmem_real_vctrs.begin();

      while (it != distmem_real_vctrs.end())
      {
        (*it)->resize(n_part);
        auto distance = std::distance(distmem_real_vctrs.begin(), it);
        thrust::copy( in_real_bfr.begin() + distance * n_copied, in_real_bfr.begin() + (distance+1) * n_copied, (*it)->begin() + n_part_old);
        it++;
      }
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::flag_lft()
    {
      thrust::copy(
        thrust::make_constant_iterator<n_t>(0),
        thrust::make_constant_iterator<n_t>(0) + lft_count,
        thrust::make_permutation_iterator(n.begin(), lft_id.begin())
      );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::flag_rgt()
    {
      thrust::copy(
        thrust::make_constant_iterator<n_t>(0),
        thrust::make_constant_iterator<n_t>(0) + rgt_count,
        thrust::make_permutation_iterator(n.begin(), rgt_id.begin())
      );
    }
  };
};
