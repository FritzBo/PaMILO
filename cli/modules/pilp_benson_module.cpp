/**
 * @file pilp_benson_module.cpp
 * @author Fritz Bökler and Mirko H. Wagner
 * @brief
 * @date 28.05.2020
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#include "pilp_benson_module.h"

using std::list;
using std::pair;
using std::string;

#include <tclap/CmdLine.h>

using TCLAP::ArgException;
using TCLAP::CmdLine;
using TCLAP::SwitchArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::ValueArg;

#include <pamilo/basic/point.h>
#include <pamilo/benchmarks/lp_parser.h>
#include <pamilo/pilp/pilp_dual_benson.h>

#include <pamilo/generic/benson_dual/ove_fp_v2.h>

#include <config_autogen.h>

using pamilo::GraphlessOVE;
using pamilo::LPparser;
using pamilo::PilpDualBensonSolver;
using pamilo::Point;

void PilpBensonModule::perform(int argc, char **argv)
{
    try
    {
        CmdLine cmd("Dual Benson to find the Pareto-frontier of the parametric integer linear "
                    "program problem.",
                    ' ', "0.1");

        ValueArg<string> output_name_argument(
            "o", "output", "Basename of the output files. This defaults to <instance>.", false, "",
            "output");

        ValueArg<double> epsilon_argument("e", "epsilon",
                                          "Epsilon to be used in floating point calculations.",
                                          false, 1E-7, "epsilon");

        ValueArg<double> point_epsilon_argument(
            "p", "point-epsilon",
            "Epsilon to decide if a potential new extreme point is already represented by an old "
            "one via euclidean distance. A value < 0 deactivates point pruning. Deacticvated by "
            "default",
            false, -1, "point-epsilon");

        ValueArg<double> solver_epsilon_argument(
            "s", "solver-epsilon",
            "Epsilon to be used in floating point calculations of the solver. This defaults to -1 "
            "(use default eps of solver).",
            false, -1, "solver-epsilon");

        ValueArg<double> vertex_enumerator_epsilon_argument(
            "v", "vertex-enumerator-epsilon",
            "Epsilon to be used in floating point calculations of the vertex enumerator. This "
            "defaults to -1 (use default epsilon of vertex enumerator).",
            false, -1, "vertex-enumerator-epsilon");

        UnlabeledValueArg<string> instance_name_argument("instance", "Name of the instance file.",
                                                         true, "", "instance");

        SwitchArg no_preprocessing_argument(
            "", "no-pre",
            "Don't run preprocessing. Only use this, if you know all objectives are in roughly the "
            "same range and either the lowest or the highest value in each objective is close to "
            "0.",
            false);

        ValueArg<string> print_type_argument("f", "solution-print-type",
                                             "Which output format for the solution file is to be "
                                             "used. Options are: json (default) and polyscip",
                                             false, "json", "solution-print-type");

        cmd.add(output_name_argument);
        cmd.add(epsilon_argument);
        cmd.add(point_epsilon_argument);
        cmd.add(solver_epsilon_argument);
        cmd.add(vertex_enumerator_epsilon_argument);
        cmd.add(instance_name_argument);
        cmd.add(no_preprocessing_argument);
        cmd.add(print_type_argument);

        cmd.parse(argc, argv);

        ILP ilp;

        string instance_name = instance_name_argument.getValue();
        double epsilon = epsilon_argument.getValue();
        double pEpsilon = point_epsilon_argument.getValue();
        double sEpsilon = solver_epsilon_argument.getValue();
        double veEpsilon = vertex_enumerator_epsilon_argument.getValue();
        if (veEpsilon = -1)
        {
            veEpsilon = epsilon;
        }
        string output_name = output_name_argument.getValue();
        bool no_preprocessing = no_preprocessing_argument.getValue();
        string solPrintType = print_type_argument.getValue();

        if (output_name == "")
        {
            output_name = instance_name;
        }
        ilp.solFile.open(output_name + "_sol");
        ilp.solPrintType = solPrintType;
        if (ilp.solPrintType == "json")
        {
            ilp.solFile << "{\n\t\"solutions\": [";
        }
        ilp.logFile.open(output_name + "_log");
        ilp.cplexFile.open(output_name + "_cplex");
        ilp.noPreprocessing = no_preprocessing;

        LPparser parser;

        try
        {
            parser.getILP(instance_name, ilp);

            {
                PilpDualBensonSolver<GraphlessOVE> solver(epsilon, pEpsilon, veEpsilon, sEpsilon);
                solver.Solve(ilp);

                solutions_.insert(solutions_.begin(), solver.solutions().cbegin(),
                                  solver.solutions().cend());
            }
        }
        catch (IloException &e)
        {
            std::cerr << e.getMessage();
        }
        if (ilp.solPrintType == "json")
        {
            ilp.solFile << "\n\t],\n\t\"solutionCount\" : " << solutions_.size() << "\n}\n";
        }
        ilp.solFile.close();
        ilp.logFile.close();
        ilp.cplexFile.close();
    }
    catch (ArgException &e)
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const string, const Point>> &PilpBensonModule::solutions()
{
    return solutions_;
}

string PilpBensonModule::statistics()
{
    string stats("");
    return stats;
}
