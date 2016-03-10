#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/tests/fig_a/",
    h5  = dir + "out_lgrngn";

  auto n = h5n(h5);

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    for (auto &plt : std::set<std::string>({"rl", "rr", "nc", "nr", "ef", "th", "rv", "u", "w"}))
    {
      Gnuplot gp;
      init(gp, h5 + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".svg", 1, 1, n); 

      if (plt == "rl")
      {
	// cloud water content
	//                                                         rho_w  kg2g
	auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
	gp << "set title 'cloud water mixing ratio [g/kg]'\n";
//	gp << "set cbrange [0:1.3]\n";
	plot(gp, tmp);
      }
      else if (plt == "rr")
      {
	// rain water content
	//                                                         rho_w  kg2g
	auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
	gp << "set logscale cb\n";
	gp << "set title 'rain water mixing ratio [g/kg]'\n";
//	gp << "set cbrange [1e-2:1]\n";
	plot(gp, tmp);
	gp << "unset logscale cb\n";
      }
      else if (plt == "nc")
      {
	// cloud particle concentration
	auto tmp = 1e-6 * h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
	gp << "set title 'cloud droplet spec. conc. [mg^{-1}]'\n";
//	gp << "set cbrange [0:150]\n";
	plot(gp, tmp);
      }
      else if (plt == "nr")
      {
	// rain particle concentration
	auto tmp = 1e-6 * h5load(h5, "rw_rng001_mom0", at * n["outfreq"]);
	gp << "set title 'rain drop spec. conc. [mg^{-1}]'\n";
//	gp << "set cbrange [.01:10]\n";
	gp << "set logscale cb\n";
	plot(gp, tmp);
	gp << "unset logscale cb\n";
      }
      else if (plt == "ef")
      {
	// effective radius
	auto r_eff = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) / h5load(h5, "rw_rng000_mom2", at * n["outfreq"]) * 1e6;
	gp << "set title 'cloud droplet effective radius [μm]'\n"; 
//	gp << "set cbrange [1:20]\n";
	plot(gp, r_eff);
      }
      else if (plt == "rv")
      {   
        // cloud particle concentration
        auto tmp = h5load(h5, "rv", at * n["outfreq"]);
        plot(gp, tmp);
      }   
      else if (plt == "th")
      {   
        // cloud particle concentration
        auto tmp = h5load(h5, "th", at * n["outfreq"]);
        plot(gp, tmp);
      }   
      else if (plt == "u")
      {   
        // cloud particle concentration
        auto tmp = h5load(h5, "u", at * n["outfreq"]);
        plot(gp, tmp);
      }   
      else if (plt == "w")
      {   
        // cloud particle concentration
        auto tmp = h5load(h5, "w", at * n["outfreq"]);
        plot(gp, tmp);
      }   

      else assert(false);
    } // var loop
  } // time loop
} // main
