/**
 * @file pilp_benson_module.cpp
 * @author Fritz Bökler and Mirko H. Wagner and Levin Nemesch
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

#ifdef USE_CDD
#    include <pamilo/generic/benson_dual/ove_cdd.h>
using pamilo::OnlineVertexEnumeratorCDD;
#endif
#include <pamilo/generic/benson_dual/ove_fp_v2.h>

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

        ValueArg<double> epsilon_argument(
            "e", "epsilon", "Epsilon to be used in floating point calculations.", false,
#ifdef USE_GRB
            1E-6
#elif USE_CPLEX
            1E-7
#endif
            ,
            "epsilon");

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

        ValueArg<string> ve_argument("E", "vertex-enumeration",
                                     "Which vertex enumeration algorithm is to be used. Options "
                                     "are: cdd, default",
                                     false, "default", "vertex-enumeration");

        ValueArg<int> solver_threads_limit(
            "t", "solver-thread-limit",
            "Maximum number of threads the solver is allowed to use. This defaults to 1.", false, 1,
            "solver-thread-limit");

        UnlabeledValueArg<string> instance_name_argument("instance", "Name of the instance file.",
                                                         true, "", "instance");

        SwitchArg no_preprocessing_argument(
            "", "no-pre",
            "Don't run preprocessing. Only use this, if you know all objectives are in roughly the "
            "same range and either the lowest or the highest value in each objective is close to "
            "0.",
            false);

        SwitchArg non_convex_argument(
            "", "non-con",
            "Allows Gurobi to also attempt solving non-convex quadratic problems. By default this "
            "is off. Non-convex problems might have especially high runtime. No effect for CPLEX.",
            false);

        ValueArg<string> print_type_argument("f", "solution-print-type",
                                             "Which output format for the solution file is to be "
                                             "used. Options are: json (default) and polyscip",
                                             false, "json", "solution-print-type");

        cmd.add(output_name_argument);
        cmd.add(epsilon_argument);
        cmd.add(point_epsilon_argument);
        cmd.add(solver_epsilon_argument);
        cmd.add(ve_argument);
        cmd.add(vertex_enumerator_epsilon_argument);
        cmd.add(solver_threads_limit);
        cmd.add(instance_name_argument);
        cmd.add(non_convex_argument);
        cmd.add(no_preprocessing_argument);
        cmd.add(print_type_argument);

        cmd.parse(argc, argv);

        ILP ilp;

#ifdef USE_GRB

        ilp.env.set(GRB_IntParam_LogToConsole, 0);
        ilp.env.set(GRB_IntParam_Threads,
                    solver_threads_limit.getValue() <= 1 ? 1 : solver_threads_limit.getValue());

        if (non_convex_argument.getValue())
        {
            ilp.env.set(GRB_IntParam_NonConvex, 2);
        }

#endif

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

        string ve = ve_argument.getValue();

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

#ifdef USE_GRB
        // Gurobi appends its logs, so we have to clean up before starting:
        std::ofstream tmpCleaner;
        tmpCleaner.open(output_name + "_gurobi");
        tmpCleaner.close();
        ilp.grbFileName = output_name + "_gurobi";
#endif

        ilp.noPreprocessing = no_preprocessing;

        LPparser parser;

        try
        {
            parser.getILP(instance_name, ilp);

#ifdef USE_CPLEX
            ilp.cplex.setParam(IloCplex::Param::Threads, solver_threads_limit.getValue() <= 1
                                                             ? 1
                                                             : solver_threads_limit.getValue());

            ilp.cplex.extract(ilp.model);
#endif

#ifdef USE_CDD
            if (ve == "cdd")
            {
                PilpDualBensonSolver<OnlineVertexEnumeratorCDD> solver(epsilon, veEpsilon,
                                                                       sEpsilon);
                solver.Solve(ilp);

                solutions_.insert(solutions_.begin(), solver.solutions().cbegin(),
                                  solver.solutions().cend());
            }
            else
#endif
            {
                if (ve == "cdd")
                {
                    std::cerr << "cdd is not activated in cmake!\n";
                    exit(0);
                }
                PilpDualBensonSolver<GraphlessOVE> solver(epsilon, veEpsilon, sEpsilon);
                solver.Solve(ilp);

                solutions_.insert(solutions_.begin(), solver.solutions().cbegin(),
                                  solver.solutions().cend());
            }
        }
#ifdef USE_GRB
        catch (GRBException &e)
        {
            if ((e.getErrorCode() == GRB_ERROR_Q_NOT_PSD) ||
                (e.getErrorCode() == GRB_ERROR_QCP_EQUALITY_CONSTRAINT))
            {
                std::cerr << "Gurobi Exception:\n"
                          << e.getMessage()
                          << "\n\nSet --non-con as argument for PaMILO to enable quadratic "
                             "non-convex optimization."
                          << std::endl;
            }
            else
            {
                std::cerr << "Gurobi Exception:\n"
                          << e.getMessage() << " with error code " << e.getErrorCode() << std::endl;
            }
        }
#elif USE_CPLEX
        catch (IloCplex::Exception &e)
        {
            std::cerr << "Cplex Exception:\n"
                      << e.getMessage() << " with status " << e.getStatus() << std::endl;
        }
#endif

        if (ilp.solPrintType == "json")
        {
            ilp.solFile << "\n\t],\n\t\"solutionCount\" : " << solutions_.size() << "\n}\n";
        }
        ilp.solFile.close();
        ilp.logFile.close();
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
