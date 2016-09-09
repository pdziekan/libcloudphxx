#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <boost/tuple/tuple.hpp>

using namespace blitz;

double iscloudy(double x)
{
  return x > 20. ? 1. : 0.;
}
BZ_DECLARE_FUNCTION(iscloudy)

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting 1 argument: dir containing out_lgrngn")

  std::string
    dir = string(av[1]),
    h5  = dir + "out_lgrngn";

  auto n = h5n(h5);
  const double z_i = 795; // [m]
  const double dz = 5; // [m], ugly and what about borde cells with dz=2.5?
  const double dx = 50; // [m], ugly and what about borde cells?
  const double dy = 1; // [m], ugly and what about borde cells?
  const double D = 3.75e-6; //[1/s], ugly, large-scale horizontal wind divergence

  Gnuplot gp;
  string file = h5 + "_series.svg";
  init_prof(gp, file, 2, 3, n); 

//  blitz::Array<float, 2> rhod;
  const double rhod=1.; //placeholder for real density table
  // read density
/*
  {
    notice_macro("about to open file: " << h5)
    H5::H5File h5f(h5 + "/const.h5", H5F_ACC_RDONLY);
  
    notice_macro("about to read dataset: G")
    H5::DataSet h5d = h5f.openDataSet("G");
    H5::DataSpace h5s = h5d.getSpace();
  
    if (h5s.getSimpleExtentNdims() != 2)  
      error_macro("need 2 dimensions")
  
    hsize_t n[2];
    enum {x, z}; 
    h5s.getSimpleExtentDims(n, NULL);
  
    rhod.resize(n[x], n[z]);
  
    hsize_t 
      cnt[2] = { n[x], n[z] },  
      off[2] = { 0,    0    };  
    h5s.selectHyperslab( H5S_SELECT_SET, cnt, off);
  
    hsize_t ext[2] = { 
      hsize_t(rhod.extent(0)), 
      hsize_t(rhod.extent(1))
    };  
    h5d.read(rhod.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(2, ext), h5s);
  }
*/

  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::Array<float, 1> res_prof(n["t"]);
  blitz::Array<float, 1> res_pos(n["t"]);
  blitz::Range all = blitz::Range::all();

  for (auto &plt : std::set<std::string>({"wvarmax", "nc", "clfrac", "lwp", "er"}))
  {
    res_prof = 0;
    res_pos = 0;

    for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
    {
      res_pos(at) = at * n["outfreq"] * n["dt"] / 3600.;
      if (plt == "clfrac")
      {
	// cloud fraction (cloudy if N_c > 20/cm^3)
        auto tmp = h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
        blitz::Array<float, 2> snap(tmp);
        snap *= rhod; // b4 it was specific moment
        snap /= 1e6; // per cm^3
        snap = iscloudy(snap);
        res_prof(at) = blitz::mean(snap); 
      }
      else if (plt == "nc")
      {
	// cloud droplet (0.5um < r < 25 um) concentration
        auto tmp = h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
        blitz::Array<float, 2> snap(tmp);
        snap *= rhod; // b4 it was specific moment
        snap /= 1e6; // per cm^3
        res_prof(at) = blitz::mean(snap); 
      }
      else if (plt == "lwp")
      {
	// liquid water path
        {
          auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp); // cloud water mixing ratio [g/kg]
          snap *= rhod; // cloud water per cubic metre (should be wet density...)
          res_prof(at) = blitz::mean(snap); 
        }
        {
          auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp); // rain water mixing ratio [g/kg]
          snap *= rhod; // rain water per cubic metre (should be wet density...)
          res_prof(at) += blitz::mean(snap); 
        }
      }
      else if (plt == "er")
      {
	//entrainment rate as in the 2009 Ackerman paper
        // to store total mixingg ratio
        blitz::Array<float, 2> rtot(n["x"], n["z"]); 
        {
          auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp); // cloud water mixing ratio [g/kg]
          rtot = snap;
        }
        {
          auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp); // rain water mixing ratio [g/kg]
          rtot += snap;
        }
        {
          auto tmp = h5load(h5, "rv", at * n["outfreq"]) * 1e3;
          blitz::Array<float, 2> snap(tmp); // vapor mixing ratio [g/kg]
          rtot += snap;
        }
        blitz::Array<int, 1> j_i(n["x"]); // index of the inversion cell
        j_i = blitz::first(rtot < 8., j);
        res_prof(at) = blitz::mean(j_i);
      }
      else if (plt == "wvarmax")
      {
        // maximum variance of vertical velocity, assume w_mean=0
        auto tmp = h5load(h5, "w", at * n["outfreq"]);
        blitz::Array<float, 2> snap(tmp);
        blitz::Array<float, 1> mean(n["z"]);
        snap = snap * snap; // 2nd power
        mean = blitz::mean(snap(j,i), j); // mean variance of w in horizontal
        res_prof(at) = blitz::max(mean);
      }  
      else assert(false);
    } // time loop

    gp << "set yrange[*:*]\n";
//    gp << "set xrange[0:6]\n";

    if (plt == "clfrac")
      gp << "set title 'average cloud fraction'\n";
    else if (plt == "nc")
      gp << "set title 'average cloud drop conc [1/cm^3]'\n";
    else if (plt == "wvarmax")
      gp << "set title 'max variance of w [m^2 / s^2]'\n";
    else if (plt == "lwp")
    {
      gp << "set title 'liquid water path [g / m^2]'\n";
      res_prof *= dz * n["z"];
    }
    else if (plt == "er")
    {
      // forward difference, in cm
      blitz::Range nolast = blitz::Range(0, n["t"]-2);
      res_prof(nolast) = (res_prof(nolast+1) - res_prof(nolast)) * dz * 1e2 / (n["dt"] * n["outfreq"]) + D * res_prof(nolast) * dz * 1e2;
      res_prof(n["t"]-1) = 0.;
      gp << "set title 'entrainment rate [cm / s]'\n";
    }
    gp << "plot '-' with line\n";
    gp.send1d(boost::make_tuple(res_pos, res_prof));

//    plot(gp, res);
  } // var loop
} // main