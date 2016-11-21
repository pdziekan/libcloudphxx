#include <libcloudph++/common/detail/fake_units.hpp>
#include <cassert>

#include <iostream>

namespace si = libcloudphxx::common::detail::fake_units::si;
using libcloudphxx::common::detail::fake_units::quantity;
using libcloudphxx::common::detail::fake_units::pow;

template <typename real_t>
void test()
{
  // empty constructor
  quantity<si::pressure, real_t> p0;

  // ctor with assignment 
  quantity<si::pressure, real_t> p1 = real_t(200) * si::pascals;
  if (p1.value != real_t(200)) assert(false),throw;

  // assignment
  p0 = p1;
  if (p0.value != p1.value) assert(false),throw;

  // multiplication
  p0 = p1 * p1;
  if (p0.value != p1.value * p1.value) assert(false),throw;

  // division
  p0 = p1 / p1;
  if (p0.value != 1) assert(false),throw;
  
  // power
  p0 = pow(p1, p0);
  if (p0.value != 200) assert(false),throw;

  // addition
  p0 = p1 + p1;
  if (p0.value != 400) assert(false),throw;

  // subtraction
  p0 = p1 - p1;
  if (p0.value != 0) assert(false),throw;

  // cast
  double a = p0;
  if (a != p0.value) assert(false),throw;

  // addition of int
  p0 = p1 * (p1 + 1);
  if (p0.value != 200*201) assert(false),throw;

  // addition to int
  p0 = p1 * (1 + p1);
  if (p0.value != 200*201) assert(false),throw;

  // int over quantity
  p0 = p1 * (1 / p1);
  if (p0.value != 200*(1./200)) assert(false),throw;

  // int times quantity
  p0 = p1 * (2 * p1);
  if (p0.value != 200*(2. * 200)) assert(false),throw;

  // int minus quantity
  p0 = p1 * (2 - p1);
  if (p0.value != 200*(2. - 200)) assert(false),throw;

  // q + int
  p0 = p1 * (p1 + 1);
  if (p0.value != 200*201) assert(false),throw;

  // q / int
  p0 = p1 * (p1 / 2);
  if (p0.value != 200*(200. / 2)) assert(false),throw;

  // q * int
  p0 = p1 * (p1 * 2);
  if (p0.value != 200*(2. * 200)) assert(false),throw;

  // quantity - int
  p0 = p1 * (p1 - 2);
  if (p0.value != 200*(200-2)) assert(false),throw;
}

int main()
{ 
  test<float>();
  test<double>();
  test<long double>();
}
