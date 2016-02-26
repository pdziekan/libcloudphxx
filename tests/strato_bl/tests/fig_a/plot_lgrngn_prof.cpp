#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <boost/tuple/tuple.hpp>

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: dir containing out_lgrngn")

  std::string
    dir = string(av[1]),
    h5  = dir + "out_lgrngn";

  auto n = h5n(h5);
  const double z_i = 795; // [m]
  const double dz = 5; // [m], ugly

  // average profile between 2h and 4h (in paper its between 2h and 6h! - do longer sims)
  int first_timestep = 7200. / n["dt"] / n["outfreq"];
  int last_timestep = 14400. /  n["dt"] / n["outfreq"];

  Gnuplot gp;
  string file = h5 + ".plot/profiles.svg";
  init_prof(gp, file, 1, 2, n); 

  for (auto &plt : std::set<std::string>({"rtot", "rliq", "thliq"}))
  {
    blitz::Array<float, 2> res(n["x"], n["z"]);
    res = 0;

    for (int at = first_timestep; at <= last_timestep; ++at) // TODO: mark what time does it actually mean!
    {
      if (plt == "rliq")
      {
	// liquid water content (cloud + rain, missing droplets with r<0.5um!)
	auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
        blitz::Array<float, 2> snap(tmp);
        res += snap;
        gp << "set title 'liquid water r [g/kg] averaged over 2h-4h, w/o rw<0.5um'\n";
      }
      else if (plt == "rtot")
      {
	// total water content (vapor + cloud + rain, missing droplets with r<0.5um!)
        {
          auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res += snap; 
        }
        {
          auto tmp = h5load(h5, "rv", at * n["outfreq"]) * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res += snap;
        }
        gp << "set title 'total water r [g/kg] averaged over 2h-4h, w/o rw<0.5um'\n";
      }
      else assert(false);
    } // time loop
    res /= last_timestep - first_timestep + 1;
    
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::Array<float, 1> res_prof(n["z"]);
    blitz::Array<float, 1> res_pos(n["z"]);
    res_pos = i * dz / z_i; 
    res_prof = blitz::mean(res(j,i), j); // average in horizontal
//    double test = blitz::mean(res(blitz::Range::all, 0));
//    std::cout << res << std::endl;
//    res_prof /= n["x"];
//    std::cout << res_prof << std::endl;
  //  std::cout << test << std::endl;
//std::cout << blitz::mean(res(blitz::Range::all(), 0)) << std::endl;
//    gp << "set xrange [0:" << tmp.extent(0)-1 << "]\n";
//    gp << "set yrange [0:" << tmp.extent(1)-1 << "]\n";
    gp << "plot '-' with line\n";
    gp.send1d(boost::make_tuple(res_prof, res_pos));

//    plot(gp, res);
  } // var loop
} // main
