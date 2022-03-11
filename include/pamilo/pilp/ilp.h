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

#include <ilcplex/ilocplex.h>

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
    IloEnv env;
    IloModel model;
    IloCplex cplex;
    IloObjective obj;
    IloObjective multiObj;
    IloNumVarArray vars;
    IloRangeArray cons;

    std::vector<double> relScale;
    std::vector<double> offset;

    std::string filename;
    std::ofstream solFile;
    std::string solPrintType;
    std::ofstream logFile;
    std::ofstream cplexFile;
    bool noPreprocessing;

    int startTime;

    int dimension;

    ILP()
        : model(env)
        , cplex(env)
        , vars(env)
        , cons(env)
        , dimension(-1)
        , filename("pilpInstance")
        , solFile("")
        , solPrintType("json")
        , logFile("")
        , cplexFile("")
        , noPreprocessing(false)
        , startTime(-1)
    {
    }

    ~ILP()
    {
        env.end();
    }
};
