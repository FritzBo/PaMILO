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
#include "pilp_benson_args.hpp"

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

#include <pamilo/pilp/ilp_interface.hpp>
#include <pamilo/pilp/grb_interface.hpp>

void PilpBensonModule::perform(int argc, char **argv)
{
    try
    {
        PilpBensonArgs args(argc, argv);

        ILP ilp;

#ifdef USE_GRB

        ilp.env.set(GRB_IntParam_LogToConsole, 0);
        ilp.env.set(GRB_IntParam_Threads,
                    args.solver_threads_limit.getValue() <= 1 ? 1 : args.solver_threads_limit.getValue());

        if (args.non_convex.getValue())
        {
            ilp.env.set(GRB_IntParam_NonConvex, 2);
        }

        ilp.env.set(GRB_IntParam_Presolve, args.grb_presovle.getValue());
        ilp.env.set(GRB_IntParam_Method, args.grb_lp.getValue());
        ilp.env.set(GRB_IntParam_ScaleFlag, args.grb_scale.getValue());

#endif

        string instance_name = args.instance_name.getValue();
        double epsilon = args.epsilon.getValue();
        double pEpsilon = args.point_epsilon.getValue();
        double sEpsilon = args.solver_epsilon.getValue();
        double veEpsilon = args.vertex_enumerator_epsilon.getValue();
        if (veEpsilon = -1)
        {
            veEpsilon = epsilon;
        }
        string output_name = args.output_name.getValue();
        bool no_preprocessing = args.no_preprocessing.getValue();
        string solPrintType = args.print_type.getValue();

        string ve = args.ve.getValue();

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
            pamilo::IlpInterface<pamilo::GRBInterface>test(args);
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
                          << "\n\nSet --non-con as arg for PaMILO to enable quadratic "
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
