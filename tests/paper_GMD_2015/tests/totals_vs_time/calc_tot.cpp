#include "../common.hpp"
#include "bins.hpp"
#include "hdf5.hpp"
#include <fstream>

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/tests/fig_a/",
    h5  = dir + "out_lgrngn";

  auto n = h5n(h5);
  std::ofstream of(string(av[1]) + "/tests/totals_vs_time/totals_vs_time.dat");
  of<<"#time\tmean cloud water content\tmean rain water content\ttotal no of SDs\tmean aerosol volume\tmean aerosol conc"<<std::endl;

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    // rain water content
    auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3. * 3.14 * 1e3 * 1e3;
    auto mean_rain_water = mean(tmp);

    // total sd number
    auto tmp2 = h5load(h5, "sd_conc", at * n["outfreq"]);
    auto sd_total = sum(tmp2);

    // cloud water content
    auto tmp3 = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3. * 3.14 * 1e3 * 1e3;
    auto mean_cloud_water = mean(tmp3);

    // mean aerosol volume per kg...
    auto tmp4 = h5load(h5, "rd_rng000_mom3", at * n["outfreq"]) * 4./3. * 3.14;
    auto mean_aero_vol = mean(tmp4);

    // mean numebr per kg..
    auto tmp5 = h5load(h5, "rd_rng000_mom0", at * n["outfreq"]);
    auto mean_aero_conc = mean(tmp5);

    of << at * n["outfreq"] 
    << "\t\t" << blitz::setw(8) << mean_cloud_water 
    << "\t\t" << blitz::setw(8) << mean_rain_water 
    << "\t\t" << blitz::setw(8) << sd_total 
    << "\t\t" << blitz::setw(8) << mean_aero_vol 
    << "\t\t" << blitz::setw(8) << mean_aero_conc
    << std::endl;
  } // time loop
  of.close();
} // main
