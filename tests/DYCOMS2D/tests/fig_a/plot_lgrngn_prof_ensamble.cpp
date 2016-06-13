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

  std::string prof_file="/out_lgrngn_profiles.dat";
  std::set<std::string> profs({"00rtot", "rliq", "thl", "wvar", "w3rd", "prflux", "clfrac", "N_c"});
  std::vector<Array<double, 1>> sums(profs.size());

  Gnuplot gp;
  std::string file = argv[1] + std::string("/out_lgrngn_mean_profiles.svg");
  auto n = h5n(argv[1] + std::string("out_lgrngn"));
  init_prof(gp, file, 3, 3, n);
  blitz::Array<float, 1> res_pos(n["z"]);
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
  res_pos = fi * n["dz"] / mean_z_i;

  int i=0;
  for (auto &plt : profs) 
  {
    if (plt == "rliq")
      gp << "set title 'liquid water r [g/kg] averaged over 2h-6h, w/o rw<0.5um'\n";
    else if (plt == "00rtot")
      gp << "set title 'total water r [g/kg] averaged over 2h-6h, w/o rw<0.5um'\n";
    else if (plt == "N_c")
      gp << "set title 'cloud droplets ( 0.5um < r < 25um) concentration [1/cm^3]'\n";
    else if (plt == "thl")
      gp << "set title 'liquid potential temp [K]'\n";
    else if (plt == "clfrac")
      gp << "set title 'cloud fraction'\n";
    else if (plt == "prflux")
      gp << "set title 'precipitation flux [W/m^2]'\n";
    else if (plt == "wvar")
      gp << "set title 'variance of w [m^2 / s^2]'\n";
    else if (plt == "w3rd")
      gp << "set title '3rd mom of w [m^3 / s^3]'\n";
    sums.at(i) /= ctr;
    gp << "plot '-' with line\n";
    gp.send1d(boost::make_tuple(sums.at(i), res_pos));
    ++i;
  }



}
