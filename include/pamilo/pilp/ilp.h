/**
 * @file ilp.h
 * @author Fritz Bökler and Mirko H. Wagner and Levin Nemesch
 * @date 20.06.2020
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#pragma once

#ifdef USE_GRB
#    include <gurobi_c++.h>
#elif USE_CPLEX
#    include <ilcplex/ilocplex.h>
#endif

#include <fstream>
#include <memory>

#include <vector>

/**
 * @brief Class to store all information about an ILP problem in. This includes the
 * solver environment, model, file information and several more
 *
 */
class ILP
{
public:
#ifdef USE_GRB
    GRBEnv env;
    std::unique_ptr<GRBModel> model;
    std::unique_ptr<GRBVar[]> vars;
    std::string grbFileName;
#elif USE_CPLEX
    IloEnv env;
    IloModel model;
    IloCplex cplex;
    IloObjective obj;
    IloObjective multiObj;
    IloNumVarArray vars;
    IloRangeArray cons;
    std::ofstream cplexFile;
#endif
    int sense_og;
    int n_vars;

    std::vector<double> relScale;
    std::vector<double> offset;

    std::string filename;
    std::ofstream solFile;
    std::string solPrintType;
    std::ofstream logFile;
    bool noPreprocessing;
    bool oneprepro;

    int startTime;

    int dimension;

    ILP()
        :
#ifdef USE_GRB
        model(nullptr)
        , vars(nullptr)
        , grbFileName("")
#elif USE_CPLEX
        model(env)
        , cplex(env)
        , vars(env)
        , cons(env)
        , cplexFile("")
#endif
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
