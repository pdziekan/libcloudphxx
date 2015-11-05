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

      // adjust opts_int for a distributed memory system
      // returns n_x_bfr
      template <class real_t>
      int distmem_opts(opts_init_t<real_t> &opts_init, const int &rank, const int &size)
      {
        int n_x_bfr = rank * get_dev_nx(opts_init, 0, size);

        opts_init.nx = detail::get_dev_nx(opts_init, rank, size);
 
        if(rank != 0)      opts_init.x0 = 0.;  // TODO: what if x0 greater than domain of first device?
        if(rank != size-1) opts_init.x1 = opts_init.nx * opts_init.dx;
        else               opts_init.x1 = opts_init.x1 - n_x_bfr * opts_init.dx;

        opts_init.n_sd_max = opts_init.n_sd_max / size + 1;

        return n_x_bfr;
      }

      template <class real_t>
      struct distmem
      {
        // host vector of the size of the total number of cells
        // only used for output by process rank==0
        thrust::host_vector<real_t> outbuf;  

        // external n_x_bfr
        int n_x_bfr;

        // a copy of the "global" opts_init
        opts_init_t<real_t> opts_init;

        int rank()
        {
#if defined(USE_MPI)
          int rank;
          MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
          return rank;
#else
          return 0;
#endif
        }

        int size()
        {
#if defined(USE_MPI)
          int size;
          MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &size));
          return size;
#else
          return 0;
#endif
        }

        bool dist_mem()
        {
#if defined(USE_MPI)
          return true;
#else
          return false;
#endif
        }

        int n_x_bfr()
        {
          return n_x_bfr + rank() * get_dev_nx(opts_init, 0, size());
        }

        opts_init_t<real_t> lcl_opts_init()
        {
          opts_init_t<real_t> lcl_opts_init(opts_init);

          lcl_opts_init.nx = get_dev_nx(opts_init, rank(), size());
   
          if(rank() != 0)        lcl_opts_init.x0 = 0.;  // TODO: what if x0 greater than domain of first device?
          if(rank() != size()-1) lcl_opts_init.x1 = lcl_opts_init.nx * lcl_opts_init.dx;
          else                   lcl_opts_init.x1 = lcl_opts_init.x1 - rank() * get_dev_nx(opts_init, 0, size()) * lcl_opts_init.dx;
  
          lcl_opts_init.n_sd_max = lcl_opts_init.n_sd_max / size() + 1;
        }

        // ctor
        distmem(const opts_init_t<real_t> &opts_init, const int &rank=-1, const int &size=-1, const int n_x_bfr=0):
          opts_init(opts_init)
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
#endif

#if defined(USE_MPI)
          int rank;
          MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
          if(rank==0)
#endif
            outbuf.resize(
              detail::m1(opts_init.nx) *
              detail::m1(opts_init.ny) *
              detail::m1(opts_init.nz));

        // init mpi
        MPI_CHECK(MPI_Init(nullptr, nullptr));
        }
      };
    }
  }
}
