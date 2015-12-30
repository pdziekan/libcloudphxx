// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

// contains definitions of members of particles_t specialized for multiple GPUs
#include <omp.h>
namespace libcloudphxx
{
  namespace lgrngn
  {
    // constructor
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::particles_t(const opts_init_t<real_t> &opts_init) :
      pdistmem(new distmem(opts_init))
    {
      // make opts_init point to global opts init
      this->opts_init = &pdistmem->opts_init;

      // opts init for this process (for MPI)
      opts_init_t<real_t> lcl_opts_init = pdistmem->lcl_opts_init();

      int dev_count;
      
      if(lcl_opts_init.src_switch) throw std::runtime_error("multi_CUDA is not yet compatible with source. Use other backend or turn off opts_init.src_switch.");
      if(lcl_opts_init.chem_switch) throw std::runtime_error("multi_CUDA is not yet compatible with chemistry. Use other backend or turn off opts_init.chem_switch.");

      // multi_CUDA works only for 2D and 3D
      if(lcl_opts_init.nz == 0)
        throw std::runtime_error("multi_CUDA backend works only for 2D and 3D simulations.");

      if (!(lcl_opts_init.x1 > lcl_opts_init.x0 && lcl_opts_init.x1 <= lcl_opts_init.nx * lcl_opts_init.dx))
        throw std::runtime_error("!(x1 > x0 & x1 <= min(1,nx)*dx)");

      // get number of available devices
      gpuErrchk(cudaGetDeviceCount(&dev_count)); 
      
      // set number of devices to use
      if(lcl_opts_init.dev_count > 0)
      {
        if(dev_count < lcl_opts_init.dev_count)
          throw std::runtime_error("number of available GPUs smaller than number of GPUs defined in opts_init");
        else 
          dev_count = lcl_opts_init.dev_count;
      }
      // copy dev_count to opts_init for threads to use
      lcl_opts_init.dev_count = dev_count;
   
      // check if all GPUs support UVA
      for (int i = 0; i < dev_count; ++i)
      {
        // Get device properties
        cudaDeviceProp devProp;
        gpuErrchk(cudaGetDeviceProperties(&devProp, i));
        if(!devProp.unifiedAddressing)
          throw std::runtime_error("All GPUs have to support Unified Virtual Addressing.");
        if(devProp.computeMode != 0)
          throw std::runtime_error("All GPUs have to be in the \"shared\" compute mode.");
      }

      // allow direct memory access between nieghbouring devices
      if(dev_count>1)
      {
        for(int dev_id = 0; dev_id < dev_count; ++dev_id)
        {
          gpuErrchk(cudaSetDevice(dev_id));
          // IDs of devices to the left/right, periodic boundary in x
          const int lft_dev = dev_id > 0 ? dev_id - 1 : lcl_opts_init.dev_count - 1,
                    rgt_dev = dev_id < lcl_opts_init.dev_count-1 ? dev_id + 1 : 0;

          // if available, allow direct memory access; otherwise copy through host memory will be done
          int can_access_peer;
          gpuErrchk(cudaDeviceCanAccessPeer(&can_access_peer, dev_id, lft_dev));
          if(can_access_peer)
            {gpuErrchk(cudaDeviceEnablePeerAccess(lft_dev, 0));}
          gpuErrchk(cudaDeviceCanAccessPeer(&can_access_peer, dev_id, rgt_dev));
          if(can_access_peer && dev_count > 2)
            {gpuErrchk(cudaDeviceEnablePeerAccess(rgt_dev, 0));}
        }
      }
      
      // resize the pointer vector
      particles.reserve(dev_count);

      // assign device to each thread and create particles_t in each
      for(int dev_id = 0; dev_id < dev_count; ++dev_id)
      {
        gpuErrchk(cudaSetDevice(dev_id));
        particles.push_back(new particles_t<real_t, CUDA>(lcl_opts_init, dev_id, dev_count, pdistmem->n_x_bfr()));
      }
    }

    // initialisation 
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::init(
      const arrinfo_t<real_t> th,
      const arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> courant_1,
      const arrinfo_t<real_t> courant_2,
      const arrinfo_t<real_t> courant_3,
      const std::map<enum chem_species_t, const arrinfo_t<real_t> > ambient_chem
    )
    {
      #pragma omp parallel num_threads(pdistmem->opts_init.dev_count)
      {
        const int dev_id = omp_get_thread_num();
        gpuErrchk(cudaSetDevice(dev_id));
        particles[dev_id].init(th, rv, rhod, courant_1, courant_2, courant_3, ambient_chem);
      }
    }
  };
};
