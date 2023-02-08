/**
 * @file grb_interface.hpp
 * @author Levin Nemesch
 * @brief
 * @date 06.02.2022
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */
#pragma once

#include <modules/pilp_benson_args.hpp>
#include <string>
#include <fstream>

namespace pamilo {
template <class SolverInterface>
struct IlpInterface {
    IlpInterface(PilpBensonArgs &args)
        : solver(args)
        , filename(args.instance_name.getValue())
        , solPrintType(args.print_type.getValue())
        , noPreprocessing(args.no_preprocessing.getValue())
        , startTime(-1)
    {
        std::string output_name = args.output_name.getValue();
        if (output_name == "")
        {
            output_name = filename;
        }
        if (solPrintType == "json")
        {
            solFile.open(output_name + "_sol.json");
        }
        else
        {
            solFile.open(output_name + "_sol");
        }
        logFile.open(output_name + "_log");
    }



    ~IlpInterface()
    {
        solFile.close();
        logFile.close();
    }

    SolverInterface solver;

    std::string filename;
    std::ofstream solFile;
    std::string solPrintType;
    std::ofstream logFile;
    bool noPreprocessing;
    int startTime;
};

}  // namespace pamilo
