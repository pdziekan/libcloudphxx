/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#include <libmpdata++/concurr/boost_thread.hpp>
#include "lgrngn_solver.hpp"
#include "setup.hpp"
#include "math.h"

using namespace libmpdataxx;

int main() 
{
  // compile-time parameters
  struct ct_params_t : ct_params_default_t
  {
    using real_t = float;
    enum { n_dims = 2 };
    enum { n_eqns = 4 };
    enum { rhs_scheme = solvers::trapez };
    enum { prs_scheme = solvers::cr };
    enum { opts = opts::nug | opts::iga | opts::fct};
    struct ix { enum {
      u, w, tht, rv,
      vip_i=u, vip_j=w, vip_den=-1
    }; };
  }; 

  using ix = typename ct_params_t::ix;

  const int r0 = 200;  // 500 // 100% humidity radius
  const int r1 = 300;  // transition region radius
  int nx = 401, ny = 401, nt = 1200; // nx, ny - no of cells, nt - time in sec
  int outfreq = 60; // in sec

  // conjugate residual
  using solver_t = lgrngn_solver<ct_params_t>;

  // run-time parameters
  solver_t::rt_params_t p;

  p.dt = .5;
  p.di = p.dj = 5.; 

  p.outfreq = int(double(outfreq) / p.dt + 0.5);
  nt = int(double(nt) / p.dt + 0.5);

  p.outdir = "wyniki/out_lgrngn";
  p.outvars = {
    {ix::u,   {.name = "u",   .unit = "m/s"}}, 
    {ix::w,   {.name = "w",   .unit = "m/s"}}, 
    {ix::tht, {.name = "tht", .unit = "K"  }},
    {ix::rv, {.name = "rv", .unit = "kg/kg"  }}
  };
  
/*  p.gnuplot_view = "map";
  p.gnuplot_output = "figure_%s_%04d.svg";
  p.gnuplot_with = "lines";
  p.gnuplot_surface = false;
//  p.gnuplot_contour = true;
//  p.gnuplot_cntrparam = "levels incremental 299.95, 0.1, 300.65";
//  p.gnuplot_cbrange = "[299.95 : 300.65]";
//  p.gnuplot_cbtics = "('299.99' 299.99, '300.10' 300.1, '300.20' 300.2, '300.30' 300.3, '300.40' 300.4, '300.50' 300.5, '300.60' 300.6)";
  p.gnuplot_palette = "defined (" 
    "299.95 '#ff0000'," //         
    "299.99 '#ff0000'," // 
    "299.99 '#ffffff'," //         /\-
    "300.00 '#ffffff'," //        /  \-
    "300.00 '#ffffff'," //  -----/    \---
    "300.05 '#ffffff'," // -----/      \---___
    "300.05 '#993399'," //     /        \-     ---
    "300.20 '#00CCFF'," //    /          \-       ---
    "300.35 '#66CC00'," //   /____________\-
    "300.50 '#FC8727'," //
    "300.65 '#FFFF00') maxcolors 14";
  p.gnuplot_term = "svg";*/
  p.prs_tol = 1e-7;
  p.grid_size = {nx, ny};

  libmpdataxx::concurr::boost_thread<
    solver_t, 
    bcond::cyclic , bcond::cyclic,
    bcond::rigid , bcond::rigid
  > slv(p);

  {
    // initial condition
    blitz::firstIndex i;
    blitz::secondIndex j;

    slv.advectee(ix::tht) = setup::th_dry()(j * p.dj);
    slv.advectee(ix::rv) = 0. + where(
      // if
      pow((i - nx/2) * p.di, 2) + 
      pow(j * p.dj - 800. , 2) <= pow(r0, 2), 
      // then
      setup::prtrb_rv()(j*p.dj), 
      // else
      setup::env_rv()(j*p.dj) +
      where(
      // if
        pow((i - nx/2) * p.di , 2) + 
        pow(j * p.dj - 800. , 2) <= pow(r1, 2), 
        // then
        (setup::prtrb_rv()(j*p.dj) - setup::env_rv()(j*p.dj)) * 
          pow(cos(M_PI/ 2. *
            (sqrt(pow((i - nx/2) * p.di,2) + pow(j * p.dj - 800.,2)) - r0) / (r1 - r0) // (r - r0) / (r1 - r0)
            ), 2),
        // else
        0.
      )
    );
 //   slv.advectee(ix::rv) = setup::rv_0;

    slv.advectee(ix::u) = 0; 
    slv.advectee(ix::w) = 0; 
    // density profile
    slv.g_factor() = setup::rhod()(j * p.dj);
  }

  // integration
  slv.advance(nt); 
};
