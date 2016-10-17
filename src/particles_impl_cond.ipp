// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <thrust/execution_policy.h>
#include <thrust/logical.h>

namespace libcloudphxx
{
  namespace lgrngn
  {

// --- Operator for testing nan values
struct isnan_test { 
    template<class real_t>
    BOOST_GPU_ENABLED bool operator()(const real_t a) const {
#if !defined(__NVCC__)
          using std::isnan;
          using std::isinf;
#endif
        return isnan(a) || isinf(a);
    }
};

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::cond(
      const real_t &dt,
      const real_t &RH_max
    ) {   

      // --- calc liquid water content before cond ---
      hskpng_sort(); 
      thrust_device::vector<real_t> &drv(tmp_device_real_cell);

      // calculating the 3rd wet moment before condensation
      moms_all();
      moms_calc(rw2.begin(), real_t(3./2.));

      // permute-copying the result to -dm_3
      // fill with 0s if not all cells will be updated in the following transform
      if(count_n!=n_cell)  thrust::fill(drv.begin(), drv.end(), real_t(0.));
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::negate<real_t>()
      );

      // calculating drop growth in a timestep using backward Euler 
      thrust::transform(
        rw2.begin(), rw2.end(),         // input - 1st arg (zip not as 1st arg not to write zip.end()
        thrust::make_zip_iterator(      // input - 2nd arg
          thrust::make_tuple(
            thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
            thrust::make_permutation_iterator(rv.begin(), ijk.begin()),
            thrust::make_permutation_iterator(T.begin(), ijk.begin()),
            thrust::make_permutation_iterator(p.begin(), ijk.begin()),
            thrust::make_permutation_iterator(RH.begin(), ijk.begin()),
            thrust::make_permutation_iterator(eta.begin(), ijk.begin()),
            rd3.begin(),
            kpa.begin(),
            vt.begin()
          )
        ), 
        rw2.begin(),                    // output
        detail::advance_rw2<real_t>(dt, RH_max)
      );

      // calculating the 3rd wet moment after condensation
      moms_calc(rw2.begin(), real_t(3./2.));

      // adding the third moment after condensation to dm_3
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // input - 2nd arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()
      );
      assert(thrust::none_of(thrust::seq, drv.begin(), drv.end(), isnan_test())); // run sequentially in this context (requires thrust 1.8.0)

      // update th and rv according to changes in third specific wet moment
      update_th_rv(drv);
    }
  };  
};
