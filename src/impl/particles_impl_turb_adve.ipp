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
      struct adve_euler
      {
        const real_t dt;
        adve_euler(const real_t &dt) : dt(dt) {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &x, const thrust::tuple<const real_t, const real_t> tpl) // tpl: 0 - turb_vel 1 - stretch
        {
          return x + thrust::get<0>(tpl) / thrust::get<1>(tpl) * dt;
        }
      };
    };

    // calc the SGS turbulent velocity component
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::turb_adve()
    {   
      thrust_device::vector<real_t> * vel_pos_a[] = {&x, &y, &z};
      std::vector<thrust_device::vector<real_t>*> vel_pos(&vel_pos_a[0], &vel_pos_a[0]+n_dims);

      thrust_device::vector<real_t> * vel_turbs_vctrs_a[] = {&up, &wp, &vp};
      std::vector<thrust_device::vector<real_t>*> vel_turbs_vctrs(&vel_turbs_vctrs_a[0], &vel_turbs_vctrs_a[0]+n_dims);

      namespace arg = thrust::placeholders;

      for(int i=0; i<n_dims-1; ++i)
      {
        thrust::transform(
          vel_pos[i]->begin(),
          vel_pos[i]->end(),
          vel_turbs_vctrs[i]->begin(),
          vel_pos[i]->begin(),
          arg::_1 + arg::_2 * opts_init.dt
        );
      }
      // stretching in z
      thrust::transform(
        vel_pos[n_dims-1]->begin(),
        vel_pos[n_dims-1]->end(),
        thrust::make_zip_iterator(   // input - 2nd arg
          thrust::make_tuple(
            vel_turbs_vctrs[n_dims-1]->begin(),
            thrust::make_permutation_iterator(vert_stretch_prof.begin(), k.begin())
          )
        ),
        vel_pos[n_dims-1]->begin(),
        detail::adve_euler<real_t>(opts_init.dt)
      );
    }
  };
};
