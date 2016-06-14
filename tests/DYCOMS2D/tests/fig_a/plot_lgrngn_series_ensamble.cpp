#include <blitz/array.h>
#include <fstream>
#include "../common.hpp"
#include "bins.hpp"
#include "gnuplot.hpp"
#include "hdf5.hpp"
#include <boost/tuple/tuple.hpp>


BZ_USING_NAMESPACE(blitz)

int main(int argc, char* argv[]) 
{
  Array<double, 1> snap;
  double mean_z_i = 0;
  int ctr = 0;

  std::string prof_file="/out_lgrngn_series.dat";
  std::set<std::string> profs({"wvarmax", "nc", "clfrac", "lwp", "er", "surf_precip", "mass_dry", "acc_precip"});
  std::vector<Array<double, 1>> sums(profs.size());

  Gnuplot gp;
  std::string file = argv[1] + std::string("/out_lgrngn_mean_series.svg");
  auto n = h5n(argv[1] + std::string("out_lgrngn"));
  init_prof(gp, file, 3, 3, n);
  blitz::Array<float, 1> res_pos(n["t"]);
  blitz::firstIndex fi;

  for(int i=1; i<argc; ++i)
  {
    std::string file = argv[i] + prof_file;
    std::cout << "reading in " << i << ": " << file << std::endl;
    ifstream iprof_file(file);
    if (iprof_file.bad())
    {
      cerr << "Unable to open file: " << file << endl;
      continue;
    }
    ctr++;
    int prof_ctr = 0;
    for (auto &plt : profs) 
    {
      iprof_file >> snap;
      if(i==1)
      {
        sums.at(prof_ctr).resize(snap.size());
        sums.at(prof_ctr)+=0;
      }
      sums.at(prof_ctr)+=snap;
      prof_ctr++;
    }
    double z_i;
    iprof_file >> z_i;
    mean_z_i += z_i;
  }
  mean_z_i /= ctr;
  res_pos = fi * n["outfreq"] * n["dt"] / 3600.;

  int i=0;
  for (auto &plt : profs) 
  {
    gp << "set yrange[*:*]\n";
    gp << "set xrange[*:*]\n";

    if (plt == "clfrac")
      gp << "set title 'average cloud fraction'\n";
    else if (plt == "nc")
      gp << "set title 'average cloud drop conc [1/cm^3]'\n";
    else if (plt == "wvarmax")
      gp << "set title 'max variance of w [m^2 / s^2]'\n";
    else if (plt == "surf_precip")
      gp << "set title 'surface precipitation [mm/d]'\n";
    else if (plt == "acc_precip")
      gp << "set title 'accumulated surface precipitation [mm]'\n";
    else if (plt == "mass_dry")
      gp << "set title 'total dry mass [g]'\n";
    else if (plt == "lwp")
      gp << "set title 'liquid water path [g / m^2]'\n";
    else if (plt == "er")
      gp << "set title 'entrainment rate [cm / s]'\n";

    sums.at(i) /= ctr;

    std::cout << plt << " " << res_pos << sums.at(i) << std::endl;
    gp << "plot '-' with line\n";
    gp.send1d(boost::make_tuple(res_pos, sums.at(i)));
    ++i;
  }



}
