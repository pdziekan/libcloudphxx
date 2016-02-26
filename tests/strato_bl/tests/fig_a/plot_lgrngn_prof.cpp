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
  const double dz = 5; // [m], ugly and what about borde cells with dz=2.5?

  // average profile between 2h and 4h (in paper its between 2h and 6h! - do longer sims)
  int first_timestep = 7200. / n["dt"] / n["outfreq"];
  int last_timestep = 14400. /  n["dt"] / n["outfreq"];

  Gnuplot gp;
  string file = h5 + ".plot/profiles.svg";
  init_prof(gp, file, 2, 2, n); 

  for (auto &plt : std::set<std::string>({"rtot", "rliq", "thliq", "wvar"}))
  {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::Array<float, 2> res(n["x"], n["z"]);
    blitz::Array<float, 1> res_prof(n["z"]);
    blitz::Array<float, 1> res_pos(n["z"]);
    blitz::Range all = blitz::Range::all();
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
      else if (plt == "wvar")
      {
	// variance of vertical velocity
	auto tmp = h5load(h5, "w", at * n["outfreq"]);
        blitz::Array<float, 2> snap(tmp);
        res_prof = blitz::mean(snap(j,i), j); // mean w in horizontal at this moment
        for(int ii = 0; ii < n["x"]; ++ii)
          snap(ii,all) = snap(ii,all) - res_prof(j); // snap is now w - w_mean
        snap = snap * snap; // 2nd power
        res += snap;
        gp << "set title 'variance of w [m^2 / s^2]'\n";
      }
      else assert(false);
    } // time loop
    res /= last_timestep - first_timestep + 1;
    
    res_pos = i * dz / z_i; 
    res_prof = blitz::mean(res(j,i), j); // average in horizontal

    gp << "plot '-' with line\n";
    gp.send1d(boost::make_tuple(res_prof, res_pos));

//    plot(gp, res);
  } // var loop
} // main
