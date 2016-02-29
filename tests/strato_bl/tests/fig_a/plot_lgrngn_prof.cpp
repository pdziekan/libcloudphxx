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
  const double dx = 50; // [m], ugly and what about borde cells?
  const double dy = 1; // [m], ugly and what about borde cells?

  // average profile between 2h and 6h (in paper its between 2h and 6h! - do longer sims)
  int first_timestep = 7200. / n["dt"] / n["outfreq"];
  int last_timestep = 21600. /  n["dt"] / n["outfreq"];

  Gnuplot gp;
  string file = h5 + ".plot/profiles.svg";
  init_prof(gp, file, 3, 2, n); 

  blitz::Array<float, 2> rhod;
  // read density
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

  for (auto &plt : std::set<std::string>({"rtot", "rliq", "wvar", "w3rd", "prflux"}))
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
        gp << "set title 'liquid water r [g/kg] averaged over 2h-6h, w/o rw<0.5um'\n";
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
        gp << "set title 'total water r [g/kg] averaged over 2h-6h, w/o rw<0.5um'\n";
      }
      else if (plt == "prflux")
      {
	// precipitation flux(doesnt include vertical volicty w!)
        { 
          auto tmp = h5load(h5, "precip_rate", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          snap = snap *  4./3 * 3.14 * 1e3 // to get mass
                     / dx / dy / dz    // averaged over cell volume, TODO: make precip rate return specific moment? wouldnt need the dx and dy
                     * 2264.76e3;      // latent heat of evaporation [J/kg]
          res += snap; 
        }
	// add vertical velocity to precipitation flux (3rd mom * w)
        { 
          auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]); // this time its a specific moment
          blitz::Array<float, 2> snap(tmp);
	  auto tmp2 = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap2(tmp2);
          snap = (snap * snap2) *  4./3 * 3.14 * 1e3 // to get mass
                     * rhod           // dry air density
                     * 2264.76e3;      // latent heat of evaporation [J/kg]
          res += snap; 
        }
        // turn 3rd mom * velocity into flux in [W/m^2]
        gp << "set title 'precipitation flux [W/m^2]'\n";
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
      else if (plt == "w3rd")
      {
	// variance of vertical velocity
	auto tmp = h5load(h5, "w", at * n["outfreq"]);
        blitz::Array<float, 2> snap(tmp);
        res_prof = blitz::mean(snap(j,i), j); // mean w in horizontal at this moment
        for(int ii = 0; ii < n["x"]; ++ii)
          snap(ii,all) = snap(ii,all) - res_prof(j); // snap is now (w - w_mean)
        snap = snap * snap * snap; // 3rd power
        res += snap;
        gp << "set title '3rd mom of w [m^3 / s^3]'\n";
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
