#pragma once

#include "thrust.hpp"

#if defined(__NVCC__)
#  include <curand.h>
#  include <limits>
#else
#  include <random>
#  include <algorithm>
#endif

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
#ifdef CURAND_H_
      // cuRAND API errors
      static const char *curandGetErrorString(curandStatus_t error)
      {
          switch (error)
          {
              case CURAND_STATUS_SUCCESS:
                  return "CURAND_STATUS_SUCCESS";
      
              case CURAND_STATUS_VERSION_MISMATCH:
                  return "CURAND_STATUS_VERSION_MISMATCH";
      
              case CURAND_STATUS_NOT_INITIALIZED:
                  return "CURAND_STATUS_NOT_INITIALIZED";
      
              case CURAND_STATUS_ALLOCATION_FAILED:
                  return "CURAND_STATUS_ALLOCATION_FAILED";
      
              case CURAND_STATUS_TYPE_ERROR:
                  return "CURAND_STATUS_TYPE_ERROR";
      
              case CURAND_STATUS_OUT_OF_RANGE:
                  return "CURAND_STATUS_OUT_OF_RANGE";
      
              case CURAND_STATUS_LENGTH_NOT_MULTIPLE:
                  return "CURAND_STATUS_LENGTH_NOT_MULTIPLE";
      
              case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED:
                  return "CURAND_STATUS_DOUBLE_PRECISION_REQUIRED";
      
              case CURAND_STATUS_LAUNCH_FAILURE:
                  return "CURAND_STATUS_LAUNCH_FAILURE";
      
              case CURAND_STATUS_PREEXISTING_FAILURE:
                  return "CURAND_STATUS_PREEXISTING_FAILURE";
      
              case CURAND_STATUS_INITIALIZATION_FAILED:
                  return "CURAND_STATUS_INITIALIZATION_FAILED";
      
              case CURAND_STATUS_ARCH_MISMATCH:
                  return "CURAND_STATUS_ARCH_MISMATCH";
      
              case CURAND_STATUS_INTERNAL_ERROR:
                  return "CURAND_STATUS_INTERNAL_ERROR";
          }
      
          return "<unknown>";
      }
#endif

      template <typename real_t, int backend>
      class rng
      {
#if !defined(__NVCC__)
	// serial version using C++11's <random>
	using engine_t = std::mt19937;
        using dist_u01_t = std::uniform_real_distribution<real_t>;
        using dist_normal01_t = std::normal_distribution<real_t>;
        using dist_un_t = std::uniform_int_distribution<unsigned int>;
	engine_t engine;
	dist_u01_t dist_u01;
	dist_normal01_t dist_normal01;
	dist_un_t dist_un;

	struct fnctr_u01
	{
          engine_t &engine;
          dist_u01_t &dist_u01;
	  real_t operator()() { return dist_u01(engine); }
	};

	struct fnctr_normal01
	{
          engine_t &engine;
          dist_normal01_t &dist_normal01;
	  real_t operator()() { return dist_normal01(engine); }
	};

	struct fnctr_un
	{
          engine_t &engine;
          dist_un_t &dist_un;
	  real_t operator()() { return dist_un(engine); }
	};

	public:

        // ctor
        rng(int seed) : engine(seed), dist_u01(0,1), dist_normal01(0,1), dist_un(0, std::numeric_limits<unsigned int>::max()) {}

	void generate_n(
	  thrust_device::vector<real_t> &u01, 
	  const thrust_size_t n
	) {
          // note: generate_n copies the third argument!!!
	  std::generate_n(u01.begin(), n, fnctr_u01({engine, dist_u01})); 
	}

	void generate_normal_n(
	  thrust_device::vector<real_t> &normal01, 
	  const thrust_size_t n
	) {
          // note: generate_n copies the third argument!!!
	  std::generate_n(normal01.begin(), n, fnctr_normal01({engine, dist_normal01})); 
	}

	void generate_n(
	  thrust_device::vector<unsigned int> &un, 
	  const thrust_size_t n
	) {
          // note: generate_n copies the third argument!!!
	  std::generate_n(un.begin(), n, fnctr_un({engine, dist_un})); 
	}
#endif
      };
 
      template <typename real_t>
      class rng<real_t, CUDA>
      {
#if defined(__NVCC__)
	// CUDA parallel version using curand

	// private member fields
	curandGenerator_t gen;
	
	public:

	rng(int seed)
	{
          {
	    curandStatus_t status = curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
          std::cerr << "curandcreategenerator return status: " << detail::curandGetErrorString(status) << std::endl;
	    assert(status == CURAND_STATUS_SUCCESS /* && "curandCreateGenerator failed"*/);
	    //_unused(status);
          }
          {
	    curandStatus_t status = curandSetPseudoRandomGeneratorSeed(gen, seed);
          std::cerr << "curandsetseed return status: " << detail::curandGetErrorString(status) << std::endl;
	    assert(status == CURAND_STATUS_SUCCESS /* && "curandSetPseudoRandomGeneratorSeed failed"*/);
            //_unused(status);
	  }
        }

	~rng()
	{
	  curandStatus_t status = curandDestroyGenerator(gen); 
          std::cerr << "curanddestroygen return status: " << detail::curandGetErrorString(status) << std::endl;
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandDestroyGenerator failed"*/);
          //_unused(status);
	}

	void generate_n(
	  thrust_device::vector<float> &v, 
	  const thrust_size_t n
	)
	{
          int dev_id;
          gpuErrchk(cudaGetDevice(&dev_id));
std::cerr << "generate_n called on device " << dev_id << std::endl;
          //thrust::fill(v.begin(), v.begin()+n, float(0.5));
          //return;
          gpuErrchk(cudaDeviceSynchronize());

          float *devData;
          gpuErrchk(cudaMalloc((void **)&devData, n*sizeof(float)));
          gpuErrchk(cudaSetDevice(dev_id));
	  curandStatus_t status = curandGenerateUniform(gen, devData, n);

          std::cerr << "curandgenerateuniform return status: " << detail::curandGetErrorString(status) << std::endl;
	  //
	  //thrust_device::vector<float> devData(n); 
	  //curandStatus_t status = curandGenerateUniform(gen, thrust::raw_pointer_cast(devData.data()), n);

	  //curandStatus_t status = curandGenerateUniform(gen, thrust::raw_pointer_cast(v.data()), n);
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
//          //_unused(status);
          gpuErrchk(cudaDeviceSynchronize());

//          thrust::copy(devData.begin(), devData.end(), v.begin());
//
std::cerr << "before copy in generate_n" << std::endl;
          gpuErrchk(cudaMemcpy(thrust::raw_pointer_cast(v.data()), devData, n * sizeof(float), cudaMemcpyDeviceToDevice));
std::cerr << "after copy in generate_n" << std::endl;
          gpuErrchk(cudaFree(devData));
          gpuErrchk(cudaDeviceSynchronize());
	}

	void generate_n(
	  thrust_device::vector<double> &v, 
	  const thrust_size_t n
	)
	{
	  curandStatus_t status = curandGenerateUniformDouble(gen, thrust::raw_pointer_cast(v.data()), n);
          std::cerr << "curandgenerateuniformdouble return status: " << detail::curandGetErrorString(status) << std::endl;
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
          //_unused(status);
          gpuErrchk(cudaDeviceSynchronize());
	}

	void generate_normal_n(
	  thrust_device::vector<float> &v, 
	  const thrust_size_t n
	)
	{
          int dev_id;
          gpuErrchk(cudaGetDevice(&dev_id));
std::cerr << "generate_normal_n called on device " << dev_id << std::endl;
          gpuErrchk(cudaDeviceSynchronize());

          float *devData;
          gpuErrchk(cudaMalloc((void **)&devData, n*sizeof(float)));
          gpuErrchk(cudaSetDevice(dev_id));
	  curandStatus_t status = curandGenerateNormal(gen, devData, n, float(0), float(1));

	  //curandStatus_t status = curandGenerateNormal(gen, thrust::raw_pointer_cast(v.data()), n, float(0), float(1));
          std::cerr << "curandgeneratenormal return status: " << detail::curandGetErrorString(status) << std::endl;
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
          //_unused(status);
std::cerr << "before copy in generate_normal_n" << std::endl;
          gpuErrchk(cudaMemcpy(thrust::raw_pointer_cast(v.data()), devData, n * sizeof(float), cudaMemcpyDeviceToDevice));
std::cerr << "after copy in generate_normal_n" << std::endl;
          gpuErrchk(cudaFree(devData));
          gpuErrchk(cudaDeviceSynchronize());
	}

	void generate_normal_n(
	  thrust_device::vector<double> &v, 
	  const thrust_size_t n
	)
	{
	  curandStatus_t status = curandGenerateNormalDouble(gen, thrust::raw_pointer_cast(v.data()), n, double(0), double(1));
          std::cerr << "curandgeneratenormaldouble return status: " << detail::curandGetErrorString(status) << std::endl;
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
          //_unused(status);
          gpuErrchk(cudaDeviceSynchronize());
	}

	void generate_n(
	  thrust_device::vector<unsigned int> &v, 
	  const thrust_size_t n
	)
	{
	  curandStatus_t status = curandGenerate(gen, thrust::raw_pointer_cast(v.data()), n);
          std::cerr << "curandgenerate unsigned int return status: " << detail::curandGetErrorString(status) << std::endl;
	  assert(status == CURAND_STATUS_SUCCESS /* && "curandGenerateUniform failed"*/);
          //_unused(status);
          gpuErrchk(cudaDeviceSynchronize());
	}
#endif
      };
    };
  };
};
