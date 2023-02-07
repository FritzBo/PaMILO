/**
 * @file grb_interface.hpp
 * @author Levin Nemesch
 * @brief
 * @date 06.02.2022
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#include <pamilo/pilp/grb_interface.hpp>

#include <algorithm>
#include <fstream>

namespace pamilo {

GRBEnv GRBInterface::init_grb_env(PilpBensonArgs &args)
{
    std::string output_name = args.output_name.getValue();
    if (output_name == "")
    {
        output_name = args.instance_name.getValue();
    }

    // Gurobi appends its logs, so we have to clean up before starting:
    std::ofstream tmpCleaner;
    tmpCleaner.open(output_name + "_gurobi");
    tmpCleaner.close();

    GRBEnv env;

    env.set(GRB_StringParam_LogFile, output_name + "_gurobi");
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_Threads,
            args.solver_threads_limit.getValue() <= 1 ? 1 : args.solver_threads_limit.getValue());

    if (args.non_convex.getValue())
    {
        env.set(GRB_IntParam_NonConvex, 2);
    }

    env.set(GRB_IntParam_Presolve, args.grb_presovle.getValue());
    env.set(GRB_IntParam_Method, args.grb_lp.getValue());
    env.set(GRB_IntParam_ScaleFlag, args.grb_scale.getValue());

    if (args.solver_epsilon.getValue() > 0)
    {
        env.set(GRB_DoubleParam_FeasibilityTol, args.solver_epsilon.getValue());
        env.set(GRB_DoubleParam_OptimalityTol, args.solver_epsilon.getValue());
    }
    else
    {
        env.set(GRB_DoubleParam_FeasibilityTol, 1E-6);
        env.set(GRB_DoubleParam_OptimalityTol, 1E-6);
    }

    return env;
}

GRBInterface::GRBInterface(PilpBensonArgs &args)
    : m_lp(init_grb_env(args), args.instance_name.getValue())
    , m_lex_lp(m_lp)
    , m_og_objectives_lp(m_lp.get(GRB_IntAttr_NumObj))
    , m_scaled_objectives_lp(m_lp.get(GRB_IntAttr_NumObj))
    , m_og_objectives_lex_lp(m_lp.get(GRB_IntAttr_NumObj))
    , m_scaled_objectives_lex_lp(m_lp.get(GRB_IntAttr_NumObj))
    , m_scale(1.0, m_lp.get(GRB_IntAttr_NumObj))
    , m_offset(0.0, m_lp.get(GRB_IntAttr_NumObj))
    , m_vars_last(m_lp.getVars())
    , m_vars_other(m_lex_lp.getVars())
    , m_m(m_lp.get(GRB_IntAttr_NumConstrs))
    , m_n(m_lp.get(GRB_IntAttr_NumVars))
    , m_d(m_lp.get(GRB_IntAttr_NumObj))
    , m_last_is_normal(true)
{
    for (int i = 0; i < m_d; i++)
    {
        m_og_objectives_lp[i] = m_lp.getObjective(i);
        m_og_objectives_lex_lp[i] = m_lex_lp.getObjective(i);
    }

    if (m_lp.get(GRB_IntAttr_ModelSense) == GRB_MINIMIZE)
    {
        m_sense_multi = Point(1.0, m_d);

        // Overwrite tolerances from lp file to 0 (which is gurobi's default value here)
        for (int i = 0; i < m_d; i++)
        {
            m_lex_lp.set(GRB_IntParam_ObjNumber, i);
            m_lex_lp.set(GRB_DoubleAttr_ObjNRelTol, 0);
            m_lex_lp.set(GRB_DoubleAttr_ObjNAbsTol, 0);

            m_scaled_objectives_lp[i] = m_og_objectives_lp[i];
            m_scaled_objectives_lex_lp[i] = m_og_objectives_lex_lp[i];
        }
    }
    else
    {
        m_sense_multi = Point(-1.0, m_d);
        m_lp.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        m_lex_lp.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        for (int i = 0; i <= m_d; i++)
        {
            m_scaled_objectives_lp[i] = -1.0 * m_og_objectives_lp[i];
            m_scaled_objectives_lex_lp[i] = -1.0 * m_og_objectives_lex_lp[i];
            m_lex_lp.set(GRB_IntParam_ObjNumber, i);
            m_lex_lp.setObjectiveN(m_scaled_objectives_lex_lp[i], i, 0, 1.0, 0.0, 0.0,
                                   m_lex_lp.get(GRB_StringAttr_ObjNName));
        }
    }

    m_offset = Point(0.0, m_d);

    m_lp.set(GRB_IntAttr_NumObj, 0);

    m_lp.update();
    m_lex_lp.update();
}

void GRBInterface::modify_objectives(const Point &scale, const Point &offset)
{
    assert(scale.dimension() == m_d && offset.dimension() == m_d);

    m_scale = scale;
    m_offset = offset;

    for (int i = 0; i < m_d; i++)
    {
        m_scaled_objectives_lp[i] =
            (m_sense_multi[i] * m_og_objectives_lp[i] + m_offset[i]) * m_scale[i];
        m_scaled_objectives_lex_lp[i] =
            (m_sense_multi[i] * m_og_objectives_lex_lp[i] + m_offset[i]) * m_scale[i];

        m_lex_lp.set(GRB_IntParam_ObjNumber, i);
        m_lex_lp.setObjectiveN(m_scaled_objectives_lex_lp[i], i, 0, 1.0, 0.0, 0.0,
                               m_lex_lp.get(GRB_StringAttr_ObjNName));
    }
}

std::pair<SolverStatus, double> GRBInterface::wss(const Point &weighting)
{
    assert(weighting.dimension() == m_d);
    GRBLinExpr scalarized_wss;
    for (int i = 0; i < m_d; i++)
    {
        scalarized_wss += m_scaled_objectives_lp[i] * weighting[i];
    }
    m_lp.setObjective(scalarized_wss, GRB_MINIMIZE);
    m_lp.update();

    m_lp.optimize();

    switch (m_lp.get(GRB_IntAttr_Status))
    {
        case GRB_OPTIMAL:
            if (!m_last_is_normal)
            {
                m_last_is_normal = true;
                std::swap(m_vars_last, m_vars_other);
            }
            return std::pair<SolverStatus, double>{SolverStatus::Success,
                                                   scalarized_wss.getValue()};
        case GRB_INFEASIBLE:
            return std::pair<SolverStatus, double>{SolverStatus::Infeasible, -1};
        case GRB_UNBOUNDED:
            return std::pair<SolverStatus, double>{SolverStatus::Unbounded, -1};
        default:
            return std::pair<SolverStatus, double>{SolverStatus::UnknownStatus, -1};
    }
}

std::pair<SolverStatus, double> GRBInterface::lex_wss(const Point &weighting,
                                                      const Point &rel_weight,
                                                      const std::vector<int> &prio)
{
    assert(weighting.dimension() == m_d);

    GRBLinExpr scalarized_wss;
    for (int i = 0; i < m_d; i++)
    {
        m_lex_lp.set(GRB_IntParam_ObjNumber, i);
        m_lex_lp.set(GRB_DoubleAttr_ObjNWeight, rel_weight[i]);
        m_lex_lp.set(GRB_IntAttr_ObjNPriority, prio[i]);

        scalarized_wss += m_scaled_objectives_lex_lp[i] * weighting[i];
    }

    m_lex_lp.update();

    m_lex_lp.optimize();

    switch (m_lex_lp.get(GRB_IntAttr_Status))
    {
        case GRB_OPTIMAL:
            if (m_last_is_normal)
            {
                m_last_is_normal = false;
                std::swap(m_vars_last, m_vars_other);
            }
            return std::pair<SolverStatus, double>{SolverStatus::Success,
                                                   scalarized_wss.getValue()};
        case GRB_INFEASIBLE:
            return std::pair<SolverStatus, double>{SolverStatus::Infeasible, -1};
        case GRB_UNBOUNDED:
            return std::pair<SolverStatus, double>{SolverStatus::Unbounded, -1};
        default:
            return std::pair<SolverStatus, double>{SolverStatus::UnknownStatus, -1};
    }
}

}  // namespace pamilo
