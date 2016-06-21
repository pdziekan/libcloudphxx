#include <cstdlib> // system()
#include <set>
#include <string>
#include <sstream> // std::ostringstream

#include "../common.hpp"
#include "bins.hpp"

using std::ostringstream;
using std::set;
using std::string;

int main(int ac, char** av)
{
  if (ac != 2) error_macro("expecting one argument - CMAKE_BINARY_DIR");

  string bins_wet_str, bins_dry_str, outdir;

  {
    ostringstream tmp;
    vector<quantity<si::length>> left_edges = bins_dry();
    for (int i = 0; i < left_edges.size()-1; ++i)
      tmp << float(left_edges[i] / si::metres) << ":" << float(left_edges[i + 1] / si::metres) << "|0;";
    bins_dry_str = tmp.str();
  }

  {
    ostringstream tmp;
    vector<quantity<si::length>> left_edges = bins_wet();
    for (int i = 0; i < left_edges.size()-1; ++i)
      tmp << float(left_edges[i] / si::metres) << ":" << float(left_edges[i + 1] / si::metres) << "|0;";
    bins_wet_str = tmp.str();
  }
/*  outdir=std::getenv("STORAGE");
  string tmp=std::getenv("NODE_CONF");
  outdir+="/"+tmp+"/out_lgrngn";
  std::cout << "output directory: " << outdir << std::endl;
*/
  string opts_common = 
    "--outfreq=300 --nt=21600 --spinup=3600 --nx=129 --ny=129 --nz=301 --dt=1 --relax_th_rv=false"; // DYCOMS: 128x300 ; 600 21600 3600
  set<string> opts_micro({
//    "--micro=blk_1m --outdir=out_blk_1m",
//    "--micro=blk_2m --outdir=out_blk_2m",
    "--adv_serial=false --async=true --micro=lgrngn --outdir=out_lgrngn --backend=multi_CUDA --sd_conc=38 --sd_conc_max=50 --sstp_cond=10 --sstp_coal=1 --dev_count=0 --rcyc=false --dev_id=-1"  
    " --gccn=false"
      " --out_wet=\""
        ".5e-6:25e-6|0,2,3;" // cloud droplets
        "25e-6:1|0,3;"         // "rain"
        "0:.5e-6|0;"         // aerosols
//        + bins_wet_str + // aerosol spectrum (wet)
        "\""
      " --out_dry=\""
        "0:1|3"
//        + bins_dry_str + // aerosol spectrum (dry)
      "\""
  });

  for (auto &opts_m : opts_micro)
  {
    ostringstream cmd;
    cmd << av[1] << "/src/icicle " << opts_common << " " << opts_m;  
    notice_macro("about to call: " << cmd.str())

    if (EXIT_SUCCESS != system(cmd.str().c_str()))
      error_macro("model run failed: " << cmd.str())
  }
}
