//
//  main.cpp
//  lib_tssb
//
//  Created by Seong-Hwan Jun on 2019-05-12.
//

#include <iostream>
#include <stdio.h>
#include <string>

#include <boost/program_options.hpp>
#include <gsl/gsl_rng.h>

#include "interface.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    string config_file;

    namespace po = boost::program_options;
    po::options_description desc("Program options");
    desc.add_options()
    ("help", "Put a help message here.")
    ("config_file,c", po::value<string>(&config_file)->required(), "path to configuration file.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    cout << "Starting program..." << endl;
    Interface interface(config_file);
    interface.Run();

    return 0;
}
