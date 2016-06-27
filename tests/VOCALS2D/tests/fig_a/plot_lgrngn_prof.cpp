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
  int k_i = 0; // inversion cell

  // average profile between 2h and 6h (in paper its between 2h and 6h! - do longer sims)
  int first_timestep = 7200. / n["dt"] / n["outfreq"];
  int last_timestep = 21600. /  n["dt"] / n["outfreq"];

  const double p_1000 = 1000.;
  const double L = 2.5e6;
  const double R_d = 287.;
  const double c_p = 1004;
  const double c_pd = c_p;

  double z_i;

  Gnuplot gp;
  string file = h5 + "_profiles.svg";
  init_prof(gp, file, 3, 3, n); 

  string prof_file = h5 + "_profiles.dat";
  std::ofstream oprof_file(prof_file);

  blitz::Array<float, 2> rhod;
  rhod.resize(n["x"], n["z"]);
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

  for (auto &plt : std::set<std::string>({"00rtot", "rliq", "thl", "wvar", "w3rd", "prflux", "clfrac", "N_c"})) // rtot has to be first
  {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::Array<float, 2> res(n["x"], n["z"]);
    blitz::Array<float, 2> res_tmp(n["x"], n["z"]);
    blitz::Array<float, 2> res_tmp2(n["x"], n["z"]);
    blitz::Array<float, 1> res_prof(n["z"]);
    blitz::Array<float, 1> res_pos(n["z"]);
    blitz::Range all = blitz::Range::all();
    res = 0;

    for (int at = first_timestep; at <= last_timestep; ++at) // TODO: mark what time does it actually mean!
    {
      if (plt == "rliq")
      {
	// liquid water content (cloud + rain, missing droplets with r<0.5um!)
        {
          auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res += snap; 
        }
        {
          auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res += snap; 
        }
        gp << "set title 'liquid water r [g/kg] averaged over 2h-6h, w/o rw<0.5um'\n";
      }
      else if (plt == "00rtot")
      {
	// total water content (vapor + cloud + rain, missing droplets with r<0.5um!)
        {
          auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res_tmp = snap; 
        }
        {
          auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3 * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res_tmp += snap; 
        }
        {
          auto tmp = h5load(h5, "rv", at * n["outfreq"]) * 1e3;
          blitz::Array<float, 2> snap(tmp);
          res_tmp += snap;
        }
        res += res_tmp;
        res_prof = blitz::mean(res_tmp(j, i), j); // average in x
        // find instantaneous inversion height
        k_i +=  blitz::first((res_prof < 6.));
        gp << "set title 'total water r [g/kg] averaged over 2h-6h, w/o rw<0.5um'\n";
      }
      else if (plt == "N_c")
      {
	// cloud drops concentration [1/cm^3]
        {
          auto tmp = h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          res += snap; 
        }
        gp << "set title 'cloud droplets ( 0.5um < r < 25um) concentration [1/cm^3]'\n";
      }
      else if (plt == "thl")
      {
	// liquid potential temp [K]
        {
          {
            auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3;
            blitz::Array<float, 2> snap(tmp);
            res_tmp2 = snap; 
          }
          {
            auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]) * 4./3 * 3.14 * 1e3;
            blitz::Array<float, 2> snap(tmp);
            res_tmp2 += snap; 
          }
          // res_tmp2 is now q_l (liq water content)
          auto tmp = h5load(h5, "th", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp); // snap is theta_dry
          res_tmp = pow(snap * pow(rhod * R_d / (p_1000 * 100), R_d / c_pd), c_pd / (c_pd - R_d)); // res_tmp is now temperature; 1 bar = 100 000Pa
          snap *= (res_tmp - res_tmp2 * L / c_p) / res_tmp; 
          res += snap; 
//          res += res_tmp2;
        }
        gp << "set title 'liquid potential temp [K]'\n";
      }
      else if (plt == "clfrac")
      {
	// cloud fraction (cloudy if N_c > 20/cm^3)
        {
          auto tmp = h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap);
          res += snap; 
        }
        gp << "set title 'cloud fraction'\n";
      }
      else if (plt == "prflux")
      {
	// precipitation flux(doesnt include vertical volicty w!)
        { 
          auto tmp = h5load(h5, "precip_rate", at * n["outfreq"]);
          blitz::Array<float, 2> snap(tmp);
          snap = snap *  4./3 * 3.14 * 1e3 // to get mass
                     / n["dx"] / n["dz"]    // averaged over cell volume, TODO: make precip rate return specific moment? wouldnt need the dx and dy
                     * 2264.76e3;      // latent heat of evaporation [J/kg]
          res += snap; 
        }
	// add vertical velocity to precipitation flux (3rd mom of cloud drops * w)
/*
        { 
          auto tmp = h5load(h5, "rw_rng000_mom3", at * n["outfreq"]); // this time its a specific moment
          blitz::Array<float, 2> snap(tmp);
	  auto tmp2 = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap2(tmp2);
          snap = - (snap * snap2) *  4./3 * 3.14 * 1e3 // to get mass
                     * rhod           // dry air density
                     * 2264.76e3;      // latent heat of evaporation [J/kg]
          res += snap; 
        }
	// add vertical velocity to precipitation flux (3rd mom of rain drops * w)
        { 
          auto tmp = h5load(h5, "rw_rng001_mom3", at * n["outfreq"]); // this time its a specific moment
          blitz::Array<float, 2> snap(tmp);
	  auto tmp2 = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 2> snap2(tmp2);
          snap = - (snap * snap2) *  4./3 * 3.14 * 1e3 // to get mass
                     * rhod           // dry air density
                     * 2264.76e3;      // latent heat of evaporation [J/kg]
          res += snap; 
        }
*/
        // turn 3rd mom * velocity into flux in [W/m^2]
        gp << "set title 'precipitation flux [W/m^2]'\n";
      }
      else if (plt == "wvar")
      {
	// variance of vertical velocity, w_mean=0
	auto tmp = h5load(h5, "w", at * n["outfreq"]);
        blitz::Array<float, 2> snap(tmp);
        snap = snap * snap; // 2nd power
        res += snap;
        gp << "set title 'variance of w [m^2 / s^2]'\n";
      }
      else if (plt == "w3rd")
      {
	// 3rd mom of vertical velocity, w_mean=0
	auto tmp = h5load(h5, "w", at * n["outfreq"]);
        blitz::Array<float, 2> snap(tmp);
        snap = snap * snap * snap; // 3rd power
        res += snap;
        gp << "set title '3rd mom of w [m^3 / s^3]'\n";
      }
      else assert(false);
    } // time loop
    res /= last_timestep - first_timestep + 1;
    
    z_i = (double(k_i)-0.5) / (last_timestep - first_timestep + 1) * n["dz"];
    std::cout << "average inversion height " << z_i;
    res_pos = (i-0.5) * n["dz"] / z_i; 
    res_prof = blitz::mean(res(j, i), j); // average in x

    gp << "plot '-' with line\n";
    gp.send1d(boost::make_tuple(res_prof, res_pos));

    oprof_file << res_prof ;

//    plot(gp, res);
  } // var loop
  oprof_file << z_i << std::endl;
} // main
