#include "common.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: CMAKE_BINARY_DIR")

  std::string
    dir = string(av[1]) + "/wyniki/",
    h5  = dir + "out_lgrngn";

  auto n = h5n(h5);

  for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
  {
    for (auto &plt : std::set<std::string>({"rv", "tht", "sd_conc", "no_acti"}))
    {
      Gnuplot gp;
      init(gp, h5 + ".plot/" + plt + "/" + zeropad(at * n["outfreq"]) + ".svg", 1, 1, n); 

      if (plt == "rv")
      {
      	// cloud particle concentration
      	auto tmp = h5load(h5, "rv", at * n["outfreq"]);
      	gp << "rv'\n";
      	plot(gp, tmp);
      }
      if (plt == "tht")
      {
      	// cloud particle concentration
      	auto tmp = h5load(h5, "tht", at * n["outfreq"]);
      	gp << "tht'\n";
      	plot(gp, tmp);
      }
      if (plt == "sd_conc")
      {
      	// cloud particle concentration
      	auto tmp = h5load(h5, "sd_conc", at * n["outfreq"]);
      	gp << "sd_conc\n";
      	plot(gp, tmp);
      }
      if (plt == "no_acti")
      {
      	// cloud particle concentration
      	auto tmp = h5load(h5, "no_acti", at * n["outfreq"]);
      	gp << "no_acti'\n";
      	plot(gp, tmp);
      }
      else assert(false);
    } // var loop
  } // time loop
} // main
