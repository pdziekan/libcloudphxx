#pragma once

#if defined(USE_MPI)
 #include <mpi.h>
#endif

// Error handling macros
#define MPI_CHECK(call) \
    if((call) != MPI_SUCCESS) { \
        std::cerr << "MPI error calling \""#call"\"\n"; \
        MPI_Abort(MPI_COMM_WORLD, -1);}

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      template<class real_t>
      int get_dev_nx(const opts_init_t<real_t> &opts_init, const int &rank, const int &size)
      {
        if(rank < size-1)
          return opts_init.nx / size + .5;
        else
          return opts_init.nx - rank * int(opts_init.nx / size + .5);
      }
    }

    template <typename real_t, backend_t device>
    struct particles_t<real_t, device>::distmem
    {
      // host vector of the size of the total number of cells
      // only used for output by process rank==0
      thrust::host_vector<real_t> outbuf;  

      // external n_x_bfr
      int n_x_bfr;

      int rank, size;

      bool dist_mem;

      // a copy of the "global" opts_init
      opts_init_t<real_t> opts_init;

      int n_x_bfr()
      {
        return n_x_bfr + rank() * detail::get_dev_nx(opts_init, 0, size());
      }

      opts_init_t<real_t> lcl_opts_init()
      {
        opts_init_t<real_t> lcl_opts_init(opts_init);

        lcl_opts_init.nx = detail::get_dev_nx(opts_init, rank(), size());
 
        if(rank() != 0)        lcl_opts_init.x0 = 0.;  // TODO: what if x0 greater than domain of first device?
        if(rank() != size()-1) lcl_opts_init.x1 = lcl_opts_init.nx * lcl_opts_init.dx;
        else                   lcl_opts_init.x1 = lcl_opts_init.x1 - rank() * detail::get_dev_nx(opts_init, 0, size()) * lcl_opts_init.dx;

        lcl_opts_init.n_sd_max = lcl_opts_init.n_sd_max / size() + 1;
      }

      // ctor
      distmem(const opts_init_t<real_t> &opts_init, const int &_rank=-1, const int &_size=-1, const int _n_x_bfr=0):
        opts_init(opts_init)
      {
        // rank and size are not externally specified - no a multi_CUDA spawn
        if(_rank==-1 && _size==-1)
        {
#if !defined(USE_MPI)
          if (
          // mpich
          std::getenv("PMI_RANK") != NULL ||
          // openmpi
          std::getenv("OMPI_COMM_WORLD_RANK") != NULL ||
          // lam
          std::getenv("LAMRANK") != NULL
          ) throw std::runtime_error("mpirun environment variable detected but libcloudphxx was compiled with MPI disabled");
          // TODO: mvapich2

          // init mpi
          MPI_CHECK(MPI_Init(nullptr, nullptr));
  
          MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
          if(rank==0)
#endif
            outbuf.resize(
              detail::m1(opts_init.nx) *
              detail::m1(opts_init.ny) *
              detail::m1(opts_init.nz));
        }
      }
    };
  }
}
