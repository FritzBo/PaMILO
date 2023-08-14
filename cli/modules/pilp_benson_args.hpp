/**
 * @file pilp_bensons.hpp
 * @author Levin Nemesch
 * @brief
 * @date 28.05.2020
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */
#pragma once

#include <tclap/CmdLine.h>

using TCLAP::ArgException;
using TCLAP::CmdLine;
using TCLAP::SwitchArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::ValueArg;

/**
 * @brief Class to allow easy exchange of args between objects or functions
 *
 */
class PilpBensonArgs
{
public:
    PilpBensonArgs(int argc, char *argv[])
    {
        cmd.add(grb_presovle);
        cmd.add(grb_lp);
        cmd.add(grb_scale);

        cmd.add(solver_choice);
        cmd.add(output_name);
        cmd.add(epsilon);
        cmd.add(point_epsilon);
        cmd.add(solver_epsilon);
        cmd.add(ve);
        cmd.add(vertex_enumerator_epsilon);
        cmd.add(solver_threads_limit);
        cmd.add(instance_name);
        cmd.add(non_convex);
        cmd.add(no_preprocessing);
        cmd.add(print_type);

        cmd.parse(argc, argv);
    }

    CmdLine cmd{"Dual Benson to find the Pareto-frontier of the parametric integer linear "
                "program problem.",
                ' ', "0.1"};

    ValueArg<std::string> output_name{
        "o",   "output", "Basename of the output files. This defaults to <instance>.",
        false, "",       "output"};

    ValueArg<double> epsilon{"e",   "epsilon", "Epsilon to be used in floating point calculations.",
                             false, 1E-6,      "epsilon"};

    ValueArg<std::string> solver_choice{
        "",
        "solver",
        "Solver to use for WSS. Choices are gurobi and cplex. Default is "
#ifdef USE_GRB
        "gurobi"
#else
        "cplex"
#endif
        ,
        false,
        "gurobi",
        "string"};

    ValueArg<double> point_epsilon{
        "p",
        "point-epsilon",
        "Epsilon to decide if a potential new extreme point is already represented by an old "
        "one via euclidean distance. A value < 0 deactivates point pruning. Deacticvated by "
        "default",
        false,
        -1,
        "point-epsilon"};

    ValueArg<double> solver_epsilon{
        "s",
        "solver-epsilon",
        "Epsilon to be used in floating point calculations of the solver. This defaults to -1 "
        "(use default eps of solver).",
        false,
        -1,
        "solver-epsilon"};

    ValueArg<double> vertex_enumerator_epsilon{
        "v",
        "vertex-enumerator-epsilon",
        "Epsilon to be used in floating point calculations of the vertex enumerator. This "
        "defaults to -1 (use default epsilon of vertex enumerator).",
        false,
        -1,
        "vertex-enumerator-epsilon"};

    ValueArg<std::string> ve{"E",
                             "vertex-enumeration",
                             "Which vertex enumeration algorithm is to be used. Options "
                             "are: cdd, default",
                             false,
                             "default",
                             "vertex-enumeration"};

    ValueArg<int> solver_threads_limit{
        "t",
        "solver-thread-limit",
        "Maximum number of threads the solver is allowed to use. This defaults to 1.",
        false,
        1,
        "solver-thread-limit"};

    UnlabeledValueArg<std::string> instance_name{"instance", "Name of the instance file.", true, "",
                                                 "instance"};

    SwitchArg no_preprocessing{
        "", "no-pre",
        "Don't run preprocessing. Only use this, if you know all objectives are in roughly the "
        "same range and either the lowest or the highest value in each objective is close to "
        "0.",
        false};

    SwitchArg non_convex{
        "", "non-con",
        "Allows Gurobi to also attempt solving non-convex quadratic problems. By default this "
        "is off. Non-convex problems might have especially high runtime. No effect for CPLEX.",
        false};

    ValueArg<std::string> print_type{"f",
                                     "solution-print-type",
                                     "Which output format for the solution file is to be "
                                     "used. Options are: json (default) and polyscip",
                                     false,
                                     "json",
                                     "solution-print-type"};

    ValueArg<int> grb_presovle{
        "", "grb-presolve", "Arg for grb presolve. Default is 0 (off).", false, 1, "int"};
    ValueArg<int> grb_lp{"",    "grb-lp", "Arg for grb lp method.Default is 0 (primal simplex).",
                         false, 0,        "int"};
    ValueArg<int> grb_scale{"", "grb-scale", "Arg for grb scaling. Default is 1.", false, 1, "int"};
};