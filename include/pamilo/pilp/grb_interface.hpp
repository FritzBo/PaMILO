/**
 * @file grb_interface.cpp
 * @author Levin Nemesch
 * @brief
 * @date 06.02.2022
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */
#pragma once

#include <memory>

#include <pamilo/basic/point.h>
#include <modules/pilp_benson_args.hpp>

#include <gurobi_c++.h>

#include "interface_util.hpp"


namespace pamilo {

class GRBInterface
{
protected:
    GRBModel m_lp;
    GRBModel m_lex_lp;

    std::vector<GRBLinExpr> m_og_objectives_lp;
    std::vector<GRBLinExpr> m_scaled_objectives_lp;
    std::vector<GRBLinExpr> m_og_objectives_lex_lp;
    std::vector<GRBLinExpr> m_scaled_objectives_lex_lp;

    Point m_sense_multi;
    Point m_scale;
    Point m_offset;

    std::unique_ptr<GRBVar[]> m_vars_last;
    std::unique_ptr<GRBVar[]> m_vars_other;

    int m_m;
    int m_n;
    unsigned int m_d;
    bool m_last_is_normal;

    static GRBEnv init_grb_env(PilpBensonArgs &args);

public:
    GRBInterface(PilpBensonArgs &args);

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
        if (m_last_is_normal)
        {
            return m_scaled_objectives_lp[d].getValue();
        }
        else
        {
            return m_scaled_objectives_lex_lp[d].getValue();
        }
    }

    inline std::string var_name(int j) const
    {
        assert(j >= 0 && j < m_n);
        return m_vars_last[j].get(GRB_StringAttr_VarName);
    }

    inline double var_value(int j) const
    {
        assert(j >= 0 && j < m_n);
        return m_vars_last[j].get(GRB_DoubleAttr_X);
    }

    inline VarType var_type(int j) const
    {
        assert(j >= 0 && j < m_n);
        if (m_vars_last[j].get(GRB_CharAttr_VType) != 'C' &&
            m_vars_last[j].get(GRB_CharAttr_VType) != 'S')
        {
            return VarType::Integer;
        }
        else
        {
            return VarType::Float;
        }
    }
};

}  // namespace pamilo