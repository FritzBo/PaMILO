/**
 * @file cplex_interface.hpp
 * @author Levin Nemesch
 * @brief
 * @date 08.02.2022
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */
#pragma once

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

#include <pamilo/basic/point.h>
#include <modules/pilp_benson_args.hpp>

#include "interface_util.hpp"

namespace pamilo {

class CPLEXInterface
{
    IloEnv m_env;
    IloModel m_model;
    IloCplex m_cplex;
    IloObjective m_og_obj;
    IloObjective m_scaled_obj;
    IloObjective m_last_obj;
    IloNumVarArray m_vars;
    IloRangeArray m_cons;
    std::ofstream m_cplexFile;

    Point m_sense_multi;
    Point m_scale;
    Point m_offset;

    int m_m;
    int m_n;
    unsigned int m_d;

public:
    CPLEXInterface(PilpBensonArgs &args);

    void modify_objectives(const Point &scale, const Point &offset);

    std::pair<SolverStatus, double> wss(const Point &weighting);

    std::pair<SolverStatus, double> lex_wss(const Point &weighting, const Point &rel_weight,
                                            const std::vector<int> &prio);

    inline int m() const
    {
        return m_m;
    }
    inline int n() const
    {
        return m_n;
    }
    inline unsigned int d() const
    {
        return m_d;
    }

    inline Point rescale(const Point &point)
    {
        assert(point.dimension() == m_d);
        Point out(m_d);
        for (int i = 0; i < m_d; i++)
        {
            out[i] = ((point[i] / m_scale[i]) - m_offset[i]) * m_sense_multi[i];
        }
        return out;
    }

    inline const Point &sense_multi() const
    {
        return m_sense_multi;
    }

    inline const Point &scale() const
    {
        return m_scale;
    }

    inline const Point &offset() const
    {
        return m_offset;
    }

    inline double obj_value(int d) const
    {
        assert(d >= 0 d < m_d);
        return m_cplex.getValue(m_scaled_obj.getCriterion(d));
    }

    inline std::string var_name(int j) const
    {
        assert(j >= 0 && j < m_n);

        return m_vars[j].getName();
    }

    inline double var_value(int j) const
    {
        assert(j >= 0 && j < m_n);

        return m_cplex.getValue(m_vars[j]);
    }

    inline VarType var_type(int j) const
    {
        assert(j >= 0 && j < m_n);
        if (m_vars[j].getType() == IloNumVar::Type::Float)
        {
            return VarType::Float;
        }
        else
        {
            return VarType::Integer;
        }
    }
};

}  // namespace pamilo
