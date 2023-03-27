/**
 * @file pilp_dual_benson.h
 * @author Fritz Bökler and Mirko H. Wagner and Levin Nemesch
 * @brief
 * @date 28.05.20
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#pragma once

#include <math.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <set>
#include <typeinfo>

#include <pamilo/pilp/ilp_interface.hpp>
#include <pamilo/pilp/interface_util.hpp>

#include <pamilo/basic/abstract_solver.h>
#include <pamilo/basic/point.h>
#include <pamilo/generic/benson_dual/dual_benson_scalarizer.h>

namespace pamilo {

/**
 * @brief Struct to count the number of calls to the ilp solver and the cumulated cpu time of those
 * calls
 *
 */
struct Log {
    /**
     * @brief Cumulated cpu time of solver calls
     *
     */
    double solver_time = 0;

    /**
     * @brief Number of solver calls
     *
     */
    int nSolves = 0;
};

/**
 * @brief Objects of this class write solutions in solution files of ILP objects. Solution type
 * used is string
 *
 */
template <class Solverinterface>
class ILPSolverPrinter
{
public:
    /**
     * @brief Construct a new ILPSolverPrinter object
     *
     * @param ilp The ILP onto whose solution file is written
     */
    ILPSolverPrinter(IlpInterface<Solverinterface> &ilp)
        : ilp_(ilp)
    {
    }

    /**
     * @brief Writes on the solution file
     *
     */
    inline void operator()(std::pair<std::string, Point *>, bool = false);

private:
    IlpInterface<Solverinterface> &ilp_;
};

template <class Solverinterface>
void ILPSolverPrinter<Solverinterface>::operator()(std::pair<std::string, Point *> sol,
                                                   bool isFirst)
{
    Point &point = *(sol.second);
    Point pointScaled(ilp_.solver.rescale((*sol.second)));

    if (ilp_.solPrintType == "polyscip")
    {
        ilp_.solFile << "[ " << pointScaled << " ]" << sol.first;
    }
    else if (ilp_.solPrintType == "json")
    {
        if (!isFirst)
        {
            ilp_.solFile << ",";
        }
        ilp_.solFile << "\n\t\t{\n\t\t\t\"values\" : [";
        for (int i = 0; i < ilp_.solver.d(); i++)
        {
            if (i != 0)
            {
                ilp_.solFile << ",";
            }
            ilp_.solFile << "\n\t\t\t\t" << pointScaled[i];
        }
        ilp_.solFile << "\n\t\t\t]," << sol.first << "\n\t\t}";
    }
#ifdef TEST
    std::cout << "1 ";
#endif
    std::cout << pointScaled << std::endl;
}

template <class Solverinterface>
class ILPSolverAdaptor
{
public:
    /**
     * @brief Construct a new ILPSolverAdaptor.
     *
     * @param ilp ILP object which contains the solver used in calculations
     * @param log The log object to which the values are added
     * @param eps Epsilon value for double comparison instead of 0.0
     */
    ILPSolverAdaptor(IlpInterface<Solverinterface> &ilp, Log &log, double eps)
        : ilp_(ilp)
        , log_(log)
        , eps_(eps)
    {
        if (ilp.noPreprocessing)
        {
            ilp.logFile << "preprocessing time: 0\n";
            return;
        }

        // Preprocessing

        clock_t start = clock();
        int dim = ilp_.solver.d();

        Point weighting(dim), mini(dim), maxi(dim);
        std::vector<double> spread(dim);
        for (int i = 0; i < dim; i++)
        {
            weighting[i] = 0;
            mini[i] = std::numeric_limits<double>::max();
            maxi[i] = std::numeric_limits<double>::lowest();
        }

        // Calculate the min and max of the vertices of the upper image for every dimension
        //  min corresponds to ideal point, max to pseudo nadir
        for (int i = 0; i < dim; i++)
        {
            weighting[i] = 1;
            Point value(dim);
            std::string sol;

            double scalar_value = operator()(weighting, value, sol);

            for (int j = 0; j < dim; j++)
            {
                mini[j] = std::min(mini[j], value[j]);
                maxi[j] = std::max(maxi[j], value[j]);
            }
            weighting[i] = 0;
        }

        int zeroCnt = 0;
        for (int i = 0; i < dim; i++)
        {
            if (maxi[i] == std::numeric_limits<double>::max() ||
                mini[i] == std::numeric_limits<double>::lowest())
            {
                std::cout << "not bounded in objective " << i << std::endl;
                exit(0);
            }
            spread[i] = maxi[i] - mini[i];
            zeroCnt += (spread[i] == 0 ? 1 : 0);
        }

        double medianSpread = 0;
        if (zeroCnt < dim)
        {
            std::nth_element(spread.begin(), spread.begin() + (dim + zeroCnt) / 2, spread.end());
            medianSpread = spread[(dim + zeroCnt) / 2];
        }

        Point p_offset(ilp_.solver.d());
        Point p_relscale(ilp_.solver.d());

        for (int i = 0; i < dim; i++)
        {
            double offset = -mini[i];
            p_offset[i] = offset;
            if (maxi[i] - mini[i] > 0 && medianSpread != 0)
            {
                p_relscale[i] = exp2(round(log2(medianSpread / (maxi[i] - mini[i]))));
            }
            else
            {
                p_relscale[i] = 1;
            }
        }

        // TODO: Make final decision which normalization to use

        ilp_.solver.modify_objectives(p_relscale, p_offset);

        ilp_.logFile << "\tRel. Scale: " << p_relscale << "\n\tOffset:     " << p_offset
                     << std::endl;

        ilp_.logFile << "preprocessing time: " << (clock() - start) / (double)CLOCKS_PER_SEC
                     << std::endl;
    }

    inline double operator()(const Point &weighting, Point &value, std::string &sol);

private:
    IlpInterface<Solverinterface> &ilp_;

    Log &log_;

    double eps_;
};

template <typename OnlineVertexEnumerator, class SolverInterface>
class PilpDualBensonSolver : public AbstractSolver<std::string>
{
public:
    /**
     * @brief Construct a new Pilp Dual Benson Solver object
     *
     * @param epsilon
     * @param pEps
     * @param veEps
     * @param solverEps
     */
    PilpDualBensonSolver(double epsilon = 1E-6, double pEps = -1, double veEps = -1)
        : epsilon_(epsilon)
        , pEps_(pEps)
        , veEps_(veEps)
    {
    }

    /**
     * @brief Creates Solution to the ILP. Writes logging messages into ilp.logFile
     *
     * @tparam OnlineVertexEnumerator The vertex enumerator to be used
     */
    void Solve(IlpInterface<SolverInterface> &ilp);

private:
    double epsilon_;
    double pEps_;
    double veEps_;
};

/**
 * @brief Solves weighted sum problem of the instance.
 *
 * @param weighting weight vector for the problem
 * @param value Point in which a solution is stored
 * @param sol String to which the point in solution space and objective space is written
 * @return double value of optimal solution
 */
template <class SolverInterface>
inline double ILPSolverAdaptor<SolverInterface>::operator()(const Point &weighting, Point &value,
                                                            std::string &sol)
{
    ilp_.logFile << "weighting " << weighting << std::endl;

    bool has_zeroes = false;
    std::vector<int> prio(ilp_.solver.d());
    Point rel_weight(ilp_.solver.d());
    for (int i = 0; i < ilp_.solver.d(); i++)
    {
        if (std::abs(weighting[i]) <= eps_)
        {
            has_zeroes = true;
            prio[i] = 0;
            rel_weight[i] = 1;
        }
        else
        {
            prio[i] = 1;
            rel_weight[i] = weighting[i];
        }
    }

    clock_t start = clock();

    std::pair<SolverStatus, double> opti_return{SolverStatus::UnknownStatus, -1};
    if (has_zeroes)
    {
        opti_return = ilp_.solver.lex_wss(weighting, rel_weight, prio);
    }
    else
    {
        opti_return = ilp_.solver.wss(weighting);
    }

    log_.solver_time += (clock() - start) / double(CLOCKS_PER_SEC);
    log_.nSolves++;

    for (int i = 0; i < ilp_.solver.d(); i++)
    {
        value[i] = ilp_.solver.obj_value(i);
    }

    // auto solveStatus = ilp_.model->get(GRB_IntAttr_Status);

    if (opti_return.first != SolverStatus::Success)
    {
        ilp_.logFile << "Subproblem could not be solved. Exiting!\n";
        exit(0);
    }

    std::stringstream sol_buffer;

    if (ilp_.solPrintType == "json")
    {
        sol_buffer << "\n\t\t\t\"variables\" : {";
    }
    bool firstPrint = true;

    // Preload values for faster iterating
    auto var_values = ilp_.solver.var_values();
    for (int i = 0; i < ilp_.solver.n(); i++)
    {
        double val = var_values[i];
        if (val > eps_ || val < -eps_)
        {
            if (ilp_.solPrintType == "json")
            {
                if (!firstPrint)
                {
                    sol_buffer << ",";
                }
                else
                {
                    firstPrint = false;
                }
                sol_buffer << "\n\t\t\t\t\"";
                sol_buffer << ilp_.solver.var_name(i);
                sol_buffer << "\" : ";
            }
            else
            {
                sol_buffer << " ";
                sol_buffer << ilp_.solver.var_name(i);
                sol_buffer << "=";
            }
            if (ilp_.solver.var_type(i) == VarType::Integer)
            {
                sol_buffer << static_cast<int>(val);
            }
            else
            {
                sol_buffer << val;
            }
        }
    }

    sol_buffer << "\n";
    if (ilp_.solPrintType == "json")
    {
        sol_buffer << "\t\t\t}";
    }
    sol += sol_buffer.str();

    return opti_return.second;
}

template <typename OnlineVertexEnumerator, class SolverInterface>
inline void PilpDualBensonSolver<OnlineVertexEnumerator, SolverInterface>::Solve(
    IlpInterface<SolverInterface> &ilp)
{
    clock_t start = clock();
    ilp.startTime = start;

    std::vector<std::pair<std::string, Point *>> frontier;
    Log log;

    DualBensonScalarizer<OnlineVertexEnumerator, std::string> dual_benson_solver(
        ILPSolverAdaptor(ilp, log, epsilon_), ILPSolverPrinter(ilp), ilp.solver.d(), epsilon_,
        pEps_, veEps_);

    dual_benson_solver.Calculate_solutions(frontier);

#ifdef TEST
    for (int i = 0; i < ilp.solver.d(); i++)
    {
        std::cout << "0";
        for (int j = 0; j < ilp.solver.d(); j++)
        {
            if (i != j)
            {
                std::cout << " 0";
            }
            else
            {
                std::cout << " 1";
            }
        }
        std::cout << std::endl;
    }
#endif

    for (auto sol : frontier)
    {
        Point pointScaled = ilp.solver.rescale(*(sol.second));
        delete sol.second;
        add_solution(sol.first, pointScaled);
    }

    ilp.logFile << "time: " << (clock() - start) / (double)CLOCKS_PER_SEC << "\n";
    ilp.logFile << "vertex enumeration time: " << dual_benson_solver.vertex_enumeration_time()
                << "\n";
    ilp.logFile << "ilp-solver time: " << log.solver_time << "\n";
    ilp.logFile << "ilp-solver solves: " << log.nSolves << "\n";
    ilp.logFile << "time per solve: " << log.solver_time / log.nSolves << "\n";
}
}  // namespace pamilo
