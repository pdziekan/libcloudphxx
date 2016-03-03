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
/*
  const double dz = 5; // [m], ugly and what about borde cells with dz=2.5?
  const double dx = 50; // [m], ugly and what about borde cells?
  const double dy = 1; // [m], ugly and what about borde cells?
*/
  Gnuplot gp;
  string file = h5 + "_series.svg";
  init_prof(gp, file, 1, 3, n); 

  blitz::Array<float, 3> rhod;
  // read density
  {
    notice_macro("about to open file: " << h5)
    H5::H5File h5f(h5 + "/const.h5", H5F_ACC_RDONLY);
  
    notice_macro("about to read dataset: G")
    H5::DataSet h5d = h5f.openDataSet("G");
    H5::DataSpace h5s = h5d.getSpace();
  
    if (h5s.getSimpleExtentNdims() != 3)  
      error_macro("need 3 dimensions")
  
    hsize_t n[3];
    enum {x, y, z}; 
    h5s.getSimpleExtentDims(n, NULL);
  
    rhod.resize(n[x], n[y], n[z]);
  
    hsize_t 
      cnt[3] = { n[x], n[y], n[z] },  
      off[3] = { 0,    0, 0    };  
    h5s.selectHyperslab( H5S_SELECT_SET, cnt, off);
  
    hsize_t ext[3] = { 
      hsize_t(rhod.extent(0)), 
      hsize_t(rhod.extent(1)), 
      hsize_t(rhod.extent(2))
    };  
    h5d.read(rhod.data(), H5::PredType::NATIVE_FLOAT, H5::DataSpace(3, ext), h5s);
  }

  for (auto &plt : std::set<std::string>({"wvarmax", "nc", "clfrac"}))
  {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;
    blitz::Array<float, 1> res_prof(n["t"]);
    blitz::Array<float, 1> res_pos(n["t"]);
    blitz::Range all = blitz::Range::all();

    for (int at = 0; at < n["t"]; ++at) // TODO: mark what time does it actually mean!
    {
      res_pos(at) = at * n["outfreq"] * n["dt"] / 3600.;
      if (plt == "clfrac")
      {
        try
        {
          // cloud fraction (cloudy if N_c > 20/cm^3)
          auto tmp = h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
          blitz::Array<float, 3> snap(tmp);
          snap *= rhod; // b4 it was specific moment
          snap /= 1e6; // per cm^3
          snap = iscloudy(snap);
          res_prof(at) = blitz::mean(snap); 
        }
        catch(...){;}
      }
      else if (plt == "nc")
      {
	// cloud droplet (0.5um < r < 25 um) concentration
        try
        {
          auto tmp = h5load(h5, "rw_rng000_mom0", at * n["outfreq"]);
          blitz::Array<float, 3> snap(tmp);
          snap /= 1e6; // per cm^3
          res_prof(at) = blitz::mean(snap); 
        }
        catch(...) {;}
      }
      else if (plt == "wvarmax")
      {
        // maximum variance of vertical velocity
        try
        {
          auto tmp = h5load(h5, "w", at * n["outfreq"]);
          blitz::Array<float, 3> snap(tmp);
          std::cout << snap << std::endl;
          blitz::Array<float, 2> mean2d(n["x"], n["z"]);
          blitz::Array<float, 1> mean(n["z"]);
          mean2d = blitz::mean(snap(i,k,j), k); // mean over second dimension (y)
          mean = blitz::mean(mean2d(j, i), j); // mean over x
          for(int ii = 0; ii < n["x"]; ++ii)
            for(int jj = 0; jj < n["y"]; ++jj)
              snap(ii, jj, all) = snap(ii, jj, all) - mean(all); // snap is now w - w_mean
          snap = snap * snap; // 2nd power
          // mean variance of w in horizontal
          mean2d = blitz::mean(snap(i,k,j), k); // mean over second dimension (y)
          mean = blitz::mean(mean2d(j, i), j); // mean over x and y
          res_prof(at) = blitz::max(mean); // the max value
        }
        catch(...) {;}
      }  
      else assert(false);
    } // time loop

    gp << "set yrange[*:*]\n";
    gp << "set xrange[0:6]\n";

    if (plt == "clfrac")
      gp << "set title 'average cloud fraction'\n";
    else if (plt == "nc")
      gp << "set title 'average cloud drop conc [1/cm^3]'\n";
    else if (plt == "wvarmax")
      gp << "set title 'max variance of w [m^2 / s^2]'\n";
    gp << "plot '-' with line\n";
    gp.send1d(boost::make_tuple(res_pos, res_prof));

//    plot(gp, res);
  } // var loop
} // main
