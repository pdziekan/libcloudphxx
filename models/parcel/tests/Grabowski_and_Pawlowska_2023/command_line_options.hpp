//
// Created by agnieszka on 19.09.23.
//

#pragma once

#include <iostream>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
namespace po = boost::program_options;

auto parse_options(int arg_count, char** arg_values) {
    po::options_description opts_main("General options");
    opts_main.add_options()
            ("help", "produce a help message")
            ("aerosol", po::value<std::string>()->required(), "one of: pristine, polluted")
            ("init", po::value<std::string>()->required(), "one of: bin, random")
            ("vertical_velocity", po::value<double>()->required(), "e.g., 0.25, 1 or 4 [m/s]")
            ("dt", po::value<double>()->required(), "e.g., 1 [s]")
            ("filename", po::value<std::string>()->required(), "name of nc file")
    ;

    po::variables_map vm;
    po::store(
            po::command_line_parser(arg_count, arg_values).options(opts_main).run(),
            vm
    );

    if (arg_count == 1 || vm.count("help")) {
        std::cout << opts_main;
        exit(EXIT_SUCCESS);
    }

    // checking if all required options present
    po::notify(vm);

    return vm;
}
