#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    using common::unary_function;

    template <typename real_t>
    struct kappa_rd_insol_t {
          real_t kappa;
          real_t rd_insol;

      kappa_rd_insol_t(real_t kappa_, real_t rd_insol_)
      : kappa(kappa_), rd_insol(rd_insol_)
      {}

      bool operator<(const kappa_rd_insol_t &other) const
      {
        if (kappa != other.kappa) return kappa < other.kappa;
        return rd_insol < other.rd_insol;
      }
        };

    // initial dry sizes of aerosol
    // defined with a distribution
    // uses shared_ptr to make opts_init copyable
    template<typename real_t>
    using dry_distros_t = std::map<
      kappa_rd_insol_t<real_t>,              // (kappa, rd_insol)
      std::shared_ptr<unary_function<real_t>> // n(ln(rd)) @ STP; alternatively it's n(ln(rd)) independent of rhod if aerosol_independent_of_rhod=true
    >;

    // defined with a size-number pair
    template<typename real_t>
    using dry_sizes_t = std::map<
      kappa_rd_insol_t<real_t>, // (kappa, rd_insol)
      std::map<real_t,           // radius [m]
        std::pair<real_t, int>   // STP_concentration [1/m^3], number of SD that represent this radius kappa and concentration
      >
    >;

    // similar, but for sources of aerosols after initialization
    template<typename real_t>
    using src_dry_distros_t = std::map<
      kappa_rd_insol_t<real_t>,              // (kappa, rd_insol)
      std::tuple<std::shared_ptr<unary_function<real_t>>, int, int> // 1st: n(ln(rd)) @ STP created per second; alternatively it's n(ln(rd)) independent of rhod if aerosol_independent_of_rhod=true; 2nd: sd_conc for this distribution ; 3rd: supstp for this aerosol (interval in timesteps beween addition of these aerosols)
    >;

    // defined with a size-number pair
    template<typename real_t>
    using src_dry_sizes_t = std::map<
      kappa_rd_insol_t<real_t>, // (kappa, rd_insol)
      std::map<real_t,           // radius [m]
        std::tuple<real_t, int, int>   // STP_concentration [1/m^3] created per second, number of SD that represent this radius kappa and concentration, supstp
      >
    >;

    // uses shared_ptr to make opts_init copyable
    // TODO: allow rd_insol (currently, its assumed to be 0)
    template<typename real_t>
    using rlx_dry_distros_t = std::unordered_map<
      real_t,                // kappa
      std::tuple<
        std::shared_ptr<unary_function<real_t>>, // n(ln(rd)) @ STP; alternatively it's n(ln(rd)) independent of rhod if aerosol_independent_of_rhod=true
        std::pair<real_t, real_t>, // kappa range of CCN considered to belong to this distribution, ranges of different members of the map need to be exclusive (TODO: add a check of this)
        std::pair<real_t, real_t>  // range of altitudes at which this relaxation acts
      >
    > ;
  };
};
