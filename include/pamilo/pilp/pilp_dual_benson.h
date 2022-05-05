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

#include <pamilo/pilp/ilp.h>

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
class ILPSolverPrinter
{
public:
    /**
     * @brief Construct a new ILPSolverPrinter object
     *
     * @param ilp The ILP onto whose solution file is written
     */
    ILPSolverPrinter(ILP &ilp)
        : ilp_(ilp)
    {
    }

    /**
     * @brief Writes on the solution file
     *
     */
    inline void operator()(std::pair<std::string, Point *>, bool);

private:
    ILP &ilp_;
};

void ILPSolverPrinter::operator()(std::pair<std::string, Point *> sol, bool isFirst = false)
{
    Point &point = *(sol.second);
    Point pointScaled(ilp_.dimension);
    for (int i = 0; i < ilp_.dimension; i++)
    {
        pointScaled[i] = (point[i] / ilp_.relScale[i]) - ilp_.offset[i];
    }
    if (ilp_.sense_og == GRB_MAXIMIZE)
    {
        pointScaled *= -1;
    }
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
        for (int i = 0; i < ilp_.dimension; i++)
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

class ILPSolverAdaptor
{
public:
    /**
     * @brief Construct a new ILPSolverAdaptor.
     *
     * @param ilp ILP object which contains the solver used in calculations
     * @param log The log object to which the values are added
     * @param eps Epsilon value for double comparison instead of 0.0
     * @param solverEps Optional epsilon value used for Gurobi linearization. Ignored if -1
     */
    ILPSolverAdaptor(ILP &ilp, Log &log, double eps, double solverEps = -1)
        : ilp_(ilp)
        , log_(log)
        , eps_(eps)
        , solverEps_(solverEps)
    {
        if (solverEps != -1)
        {
            ilp.logFile << "solver epsilon set to: " << solverEps << "\n";
            ilp.model->set(GRB_DoubleParam_FeasibilityTol, solverEps_);
            ilp.model->update();
        }
        if (ilp.noPreprocessing)
        {
            ilp.logFile << "preprocessing time: 0\n";
            return;
        }

        // Preprocessing

        clock_t start = clock();
        int dim = ilp_.dimension;

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

        //todo roll back or delete
        // double medianSpread = 0;
        // if (zeroCnt < dim)
        // {
        //     std::nth_element(spread.begin(), spread.begin() + (dim + zeroCnt) / 2, spread.end());
        //     medianSpread = spread[(dim + zeroCnt) / 2];
        // }

        // for (int i = 0; i < dim; i++)
        // {
        //     double offset = -mini[i];
        //     ilp_.offset[i] = offset;
        //     if (maxi[i] - mini[i] > 0 && medianSpread != 0)
        //     {
        //         ilp_.relScale[i] = exp2(round(log2(medianSpread / (maxi[i] - mini[i]))));
        //     }
        // }

        for (int i = 0; i < dim; i++)
        {
            double offset = -mini[i];
            ilp_.offset[i] = offset;
            if (maxi[i] - mini[i] > 0)
            {
                ilp_.relScale[i] = exp2(round(log2(1 / (maxi[i] - mini[i]))));
            }
        }

        for (int i = 0; i < ilp_.dimension; i++)
        {
            ilp_.model->setObjectiveN(
                ilp_.relScale[i] * (ilp_.model->getObjective(i) + ilp_.offset[i]), i);
        }

        ilp_.model->update();

        ilp_.logFile << "preprocessed objectives:\n";
        for (int i = 0; i < ilp_.dimension; i++)
        {
            ilp_.logFile << "\t" << ilp_.model->getObjective(i) << "\n";
        }
        ilp_.logFile << "preprocessing time: " << (clock() - start) / (double)CLOCKS_PER_SEC
                     << std::endl;
    }

    inline double operator()(const Point &weighting, Point &value, std::string &sol);

private:
    ILP &ilp_;

    Log &log_;

    double eps_;
    double solverEps_;
};

template <typename OnlineVertexEnumerator>
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
    PilpDualBensonSolver(double epsilon = 1E-7, double pEps = 1E-5, double veEps = -1, double solverEps = -1)
        : epsilon_(epsilon)
        , pEps_(pEps)
        , veEps_(veEps)
        , solverEps_(solverEps)
    {
    }

    /**
     * @brief Creates Solution to the ILP. Writes logging messages into ilp.logFile
     *
     * @tparam OnlineVertexEnumerator The vertex enumerator to be used
     */
    void Solve(ILP &ilp);

private:
    double epsilon_;
    double pEps_;
    double veEps_;
    double solverEps_;
};

/**
 * @brief Solves weighted sum problem of the instance.
 *
 * @param weighting weight vector for the problem
 * @param value Point in which a solution is stored
 * @param sol String to which the point in solution space and objective space is written
 * @return double value of optimal solution
 */
inline double ILPSolverAdaptor::operator()(const Point &weighting, Point &value, std::string &sol)
{
    ilp_.logFile << "weighting " << weighting << std::endl;
    
    for (int i = 0; i < ilp_.dimension; i++)
    {
        ilp_.model->set(GRB_IntParam_ObjNumber, i);
        if (weighting[i] >= eps_)
        {
            ilp_.model->set(GRB_DoubleAttr_ObjNWeight, weighting[i]);
            ilp_.model->set(GRB_IntAttr_ObjNPriority, 1);
        }
        else
        {
            ilp_.model->set(GRB_DoubleAttr_ObjNWeight, 1);
            ilp_.model->set(GRB_IntAttr_ObjNPriority, 0);
        }
    }

    ilp_.model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    ilp_.model->update();

    bool doRun = true;
    clock_t start = clock();

    ilp_.model->optimize();

    log_.solver_time += (clock() - start) / double(CLOCKS_PER_SEC);
    log_.nSolves++;

    for (int i = 0; i < ilp_.dimension; i++)
    {
        value[i] = ilp_.model->getObjective(i).getValue();
    }

    auto solveStatus = ilp_.model->get(GRB_IntAttr_Status);

    if (solveStatus == GRB_INF_OR_UNBD || solveStatus == GRB_INFEASIBLE ||
        solveStatus == GRB_UNBOUNDED)
    {
        ilp_.logFile << "Problem is unbounded or infeasible!\n";
        exit(0);
    }

    if (ilp_.solPrintType == "json")
    {
        sol += "\n\t\t\t\"variables\" : {";
    }
    bool firstPrint = true;
    for (int i = 0; i < ilp_.n_vars; i++)
    {
        const auto &var = ilp_.vars[i];
        float val = var.get(GRB_DoubleAttr_X);
        if (val > eps_ || val < -eps_)
        {
            if (ilp_.solPrintType == "json")
            {
                if (!firstPrint)
                {
                    sol += ",";
                }
                else
                {
                    firstPrint = false;
                }
                sol += "\n\t\t\t\t\"";
                sol += var.get(GRB_StringAttr_VarName);
                sol += "\" : ";
            }
            else
            {
                sol += " ";
                sol += var.get(GRB_StringAttr_VarName);
                sol += "=";
            }
            if ((var.get(GRB_CharAttr_VType) != 'C' && var.get(GRB_CharAttr_VType) != 'S') ||
                (val + eps_ > int(val - eps_) && val - eps_ < int(val + eps_)))
            {
                sol += std::to_string(int(val));
            }
            else
            {
                sol += std::to_string(val);
            }
        }

        // Todo: Test warm start for MIPs
        // Experimental: Enables warm start, if one solution is known
        // if (ilp_.model->get(GRB_IntAttr_IsMIP))
        // {
        //     ilp_.vars[i].set(GRB_DoubleAttr_Start, ilp_.vars[i].get(GRB_DoubleAttr_X));
        // }
    }
    sol += "\n";
    if (ilp_.solPrintType == "json")
    {
        sol += "\t\t\t}";
    }

    // todo
    // std::cout << std::setprecision(std::numeric_limits<double>::digits10) << "MaxVio: " << ilp_.model->get(GRB_DoubleAttr_MaxVio) << "\tCon: " << ilp_.model->get(GRB_DoubleAttr_ConstrVio) << "\tBound: " << ilp_.model->get(GRB_DoubleAttr_BoundVio) << "\tInt: " << ilp_.model->get(GRB_DoubleAttr_IntVio) << "\n";

    double scalar{0.0};
    
    for (int i = 0; i < ilp_.dimension; i++)
    {
        if (weighting[i] >= eps_)
        {
            scalar += weighting[i] * ilp_.model->getObjective(i).getValue();
        }
    }

    return scalar;

    //return singleObj.getValue();
}

template <typename OnlineVertexEnumerator>
inline void PilpDualBensonSolver<OnlineVertexEnumerator>::Solve(ILP &ilp)
{
    clock_t start = clock();
    ilp.startTime = start;

    std::vector<std::pair<std::string, Point *>> frontier;
    Log log;

    DualBensonScalarizer<OnlineVertexEnumerator, std::string> dual_benson_solver(
        ILPSolverAdaptor(ilp, log, epsilon_, solverEps_), ILPSolverPrinter(ilp), ilp.dimension,
        epsilon_, pEps_, veEps_);

    dual_benson_solver.Calculate_solutions(frontier);

#ifdef TEST
    for (int i = 0; i < ilp.dimension; i++)
    {
        std::cout << "0";
        for (int j = 0; j < ilp.dimension; j++)
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
        Point &point = *(sol.second);
        Point pointScaled(ilp.dimension);
        for (int i = 0; i < ilp.dimension; i++)
        {
            pointScaled[i] = (point[i] / ilp.relScale[i]) - ilp.offset[i];
        }
        delete sol.second;
        add_solution(sol.first, pointScaled);
    }

    ilp.logFile << "time: " << (clock() - start) / (double)CLOCKS_PER_SEC << "\n";
    ilp.logFile << "vertex enumeration time: " << dual_benson_solver.vertex_enumeration_time()
                << "\n";
    ilp.logFile << "gurobi time: " << log.solver_time << "\n";
    ilp.logFile << "gurobi solves: " << log.nSolves << "\n";
    ilp.logFile << "time per solve: " << log.solver_time / log.nSolves << "\n";
}
}  // namespace pamilo
