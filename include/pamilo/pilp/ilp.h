/**
 * @file ilp.h
 * @author Fritz Bökler
 * @brief Mirko H. Wagner
 * @date 20.06.2020
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#pragma once

#include <gurobi_c++.h>
#include <fstream>
#include <memory>

#include <vector>

/**
 * @brief Class to store all information about an ILP problem in. This includes the
 * CPLEX model, file information and several more
 *
 *
 */
class ILP
{
public:
    GRBEnv env;
    std::unique_ptr<GRBModel> model;
    int sense_og;
    GRBVar* vars;
    int n_vars;

    std::vector<double> relScale;
    std::vector<double> offset;

    std::string filename;
    std::ofstream solFile;
    std::string solPrintType;
    std::ofstream logFile;
    std::string grbFileName;
    bool noPreprocessing;

    int startTime;

    int dimension;

    ILP()
        : model(nullptr)
        , vars(nullptr)
        , n_vars(-1)
        , dimension(-1)
        , filename("pilpInstance")
        , solFile("")
        , solPrintType("json")
        , logFile("")
        , noPreprocessing(false)
        , startTime(-1)
    {
    }

    ~ILP()
    {
    }
};
