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
#include <pamilo/pilp/pilp_dual_benson.h>

#ifdef USE_CDD
#    include <pamilo/generic/benson_dual/ove_cdd.h>
using pamilo::OnlineVertexEnumeratorCDD;
#endif
#include <pamilo/generic/benson_dual/ove_fp_v2.h>

using pamilo::GraphlessOVE;

using pamilo::PilpDualBensonSolver;
using pamilo::Point;

#ifdef USE_CPLEX
#include <pamilo/pilp/cplex_interface.hpp>
#endif

#ifdef USE_GRB
#include <pamilo/pilp/grb_interface.hpp>
#endif

#include <pamilo/pilp/ilp_interface.hpp>

/**
 * @brief Templated helper to avoid spaghetti
 *
 * @tparam SolverInterface
 * @param args
 */
template <class SolverInterface>
inline void start_algo(PilpBensonArgs &args, list<pair<const string, const Point>> &solutions)
{
    auto ve = args.ve.getValue();
    auto epsilon = args.epsilon.getValue();
    auto veEpsilon = args.vertex_enumerator_epsilon.getValue();

    pamilo::IlpInterface<SolverInterface> ilp(args);

    if (ilp.solPrintType == "json")
    {
        ilp.solFile << "{\n\t\"solutions\": [";
    }

#ifdef USE_CDD
    if (ve == "cdd")
    {
        PilpDualBensonSolver<OnlineVertexEnumeratorCDD, SolverInterface> solver(epsilon,
                                                                                     veEpsilon);
        solver.Solve(ilp);

        solutions.insert(solutions.begin(), solver.solutions().cbegin(),
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

        PilpDualBensonSolver<GraphlessOVE, SolverInterface> solver(epsilon, veEpsilon);

        solver.Solve(ilp);

        solutions.insert(solutions.begin(), solver.solutions().cbegin(), solver.solutions().cend());

        if (ilp.solPrintType == "json")
        {
            ilp.solFile << "\n\t],\n\t\"solutionCount\" : " << solutions.size() << "\n}\n";
        }
    }
}

void PilpBensonModule::perform(int argc, char **argv)
{
    try
    {
        PilpBensonArgs args(argc, argv);

        try
        {
            if (args.solver_choice.getValue() == "gurobi")
            {
#ifndef USE_GRB
                std::cerr << "Gurobi is not enabled in CMake!" << std::endl;
                exit(-1);
#else
                start_algo<pamilo::GRBInterface>(args, solutions_);
#endif
            }
            else if (args.solver_choice.getValue() == "cplex")
            {
#ifndef USE_CPLEX
                std::cerr << "CPLEX is not enabled in CMake!" << std::endl;
                exit(-1);
#else
                start_algo<pamilo::CPLEXInterface>(args, solutions_);
#endif
            }
            else
            {
                std::cerr << "No valid solver: " << args.solver_choice.getValue() << std::endl;
                exit(-1);
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
#endif
#ifdef USE_CPLEX
        catch (IloCplex::Exception &e)
        {
            std::cerr << "Cplex Exception:\n"
                      << e.getMessage() << " with status " << e.getStatus() << std::endl;
        }
#endif
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
