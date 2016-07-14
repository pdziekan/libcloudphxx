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
  Array<double, 1> res;
  blitz::firstIndex fi;

  std::string prof_file_name="/out_lgrngn_mean_profiles.dat";
  std::set<std::string> profs({"00rtot", "rliq", "thl", "wvar", "w3rd", "prflux", "clfrac", "N_c", "act_rd", "gccn_rw", "gccn_rw_down"});

  blitz::Array<float, 1> res_pos;

  int prof_ctr = 0;
  for (auto &plt : profs) 
  {
    Gnuplot gp;
    std::string file = argv[1] + std::string("/profiles_compare") + plt + std::string(".svg");
    std::cout << "out file: " << file << std::endl;
    init_prof(gp, file, 1, 1);
    gp << "set ylabel 'height / inversion height'\n";

    if (plt == "rliq")
      gp << "set xlabel 'liquid water mixing ratio [g/kg]'\n";
    else if (plt == "00rtot")
      gp << "set xlabel 'total water mixing ratio [g/kg]'\n";
    else if (plt == "N_c")
      gp << "set xlabel 'cloud droplets concentration [1/cm^3]'\n";
    else if (plt == "thl")
      gp << "set xlabel 'liquid potential temp [K]'\n";
    else if (plt == "clfrac")
      gp << "set xlabel 'cloud fraction'\n";
    else if (plt == "prflux")
      gp << "set xlabel 'precipitation flux [W/m^2]'\n";
    else if (plt == "wvar")
      gp << "set xlabel 'variance of w [m^2 / s^2]'\n";
    else if (plt == "w3rd")
      gp << "set xlabel '3rd mom of w [m^3 / s^3]'\n";
    else if (plt == "act_rd")
    {
      gp << "set xlabel 'activated droplets mean dry radius [um]'\n";
      gp << "set xrange [0:0.15]\n";
    }
    else if (plt == "gccn_rw")
      gp << "set xlabel 'GCCN-based droplets mean wet radius [um]'\n";
    else if (plt == "gccn_rw_down")
      gp << "set xlabel 'GCCN-based droplets mean wet radius [um] in downdraught regions'\n";

    gp << "plot '-' with line t '" << argv[2] << "'";
    for(int j=3; j<argc; j+=2)
      gp << ", '-' with line t '" << argv[j+1] << "'";
    gp << "\n";

    for(int i=1; i<argc; i+=2)
    {
      std::string file = argv[i] + prof_file_name;
      std::cout << "reading in " << i << ": " << file << std::endl;
      ifstream iprof_file(file);
      if (iprof_file.bad())
      {
        cerr << "Unable to open file: " << file << endl;
        continue;
      }
      for(int k=0; k<=prof_ctr; ++k)
        iprof_file >> snap;

      res.resize(snap.shape());
      res_pos.resize(snap.shape());
      res = snap;

      for(int k=prof_ctr+1; k<profs.size(); ++k)
        iprof_file >> snap;
      double z_i;
      iprof_file >> z_i;
      res_pos = (fi-1.5) * 5. / z_i; // hardcoded dz = 5m !!
      std::cout << i << " " << plt << " res: " << res;
      std::cout << i << " " << plt << " res_pos: " << res_pos;
      gp.send1d(boost::make_tuple(res, res_pos));
    }
    prof_ctr++;
  }
}
