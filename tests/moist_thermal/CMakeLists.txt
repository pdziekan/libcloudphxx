if(APPLE)
  # needed for the XCode clang to be identified as AppleClang and not Clang
  cmake_minimum_required(VERSION 3.0)
else()
  # needed for the OpenMP test to work in C++-only project 
  # (see http://public.kitware.com/Bug/view.php?id=11910)
  cmake_minimum_required(VERSION 2.8.8)
endif()

project(moist_thermal CXX)

find_package(libmpdata++)

find_package(HDF5 COMPONENTS CXX REQUIRED QUIET)
find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)

if(NOT libmpdataxx_FOUND)
  message(FATAL_ERROR "local libmpdata++-config.cmake not found!")
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${libmpdataxx_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELEASE "")
else()
  set(CMAKE_CXX_FLAGS_DEBUG ${libmpdataxx_CXX_FLAGS_DEBUG})
endif()

# to make <libmpdata++/...> work
set(CMAKE_CXX_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CXX_FLAGS}")

# search for a local version of libcloudphxx
find_library(libcloudphxx_LIBRARIES cloudphxx_lgrngn ${CMAKE_SOURCE_DIR}/../../build/src NO_DEFAULT_PATH)
if (libcloudphxx_LIBRARIES)
  message("icicle will use libclouphxx++ from current repository")
  include_directories(PUBLIC ${CMAKE_SOURCE_DIR}/../../include)
else()
  message("libcloudphxx++ not found in current repository, icicle will use an installed version")
endif()

# if no local version found, look for an installed one
find_library(libcloudphxx_LIBRARIES cloudphxx_lgrngn)

if (NOT libcloudphxx_LIBRARIES)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()

# std::future
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <future>\n int main() { std::future<void> f; }\n" STD_FUTURE_WORKS)
if(${STD_FUTURE_WORKS})
  add_definitions(-DSTD_FUTURE_WORKS)
endif()

enable_testing()

add_executable(moist_thermal moist_thermal.cpp)

target_link_libraries(moist_thermal ${libmpdataxx_LIBRARIES})
target_include_directories(moist_thermal PUBLIC ${libmpdataxx_INCLUDE_DIRS})

find_package(Boost COMPONENTS thread iostreams system timer program_options filesystem REQUIRED)
target_link_libraries(moist_thermal ${Boost_LIBRARIES})

target_link_libraries(moist_thermal ${libcloudphxx_LIBRARIES})

add_test(moist_thermal moist_thermal)

add_executable(plot_lgrngn plot_lgrngn.cpp)
add_test(plot_lgrngn plot_lgrngn ${CMAKE_BINARY_DIR})
target_link_libraries(plot_lgrngn ${HDF5_LIBRARIES})
target_include_directories(plot_lgrngn PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(plot_lgrngn ${Boost_LIBRARIES})
target_include_directories(plot_lgrngn PUBLIC ${Boost_INCLUDE_DIRS})
