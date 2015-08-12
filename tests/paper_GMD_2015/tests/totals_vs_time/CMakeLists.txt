find_package(HDF5 COMPONENTS CXX REQUIRED QUIET)
find_package(Boost COMPONENTS iostreams system filesystem REQUIRED)

add_executable(calc_tot calc_tot.cpp)
add_test(calc_tot calc_tot ${CMAKE_BINARY_DIR})
target_link_libraries(calc_tot ${HDF5_LIBRARIES})
target_include_directories(calc_tot PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(calc_tot ${Boost_LIBRARIES})
target_include_directories(calc_tot PUBLIC ${Boost_INCLUDE_DIRS})
