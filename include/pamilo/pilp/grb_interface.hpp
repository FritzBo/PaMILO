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

#ifdef USE_GRB

#    include <fstream>
#    include <memory>

#    include <pamilo/basic/point.h>
#    include <modules/pilp_benson_args.hpp>

#    include <gurobi_c++.h>

#    include "interface_util.hpp"

namespace pamilo {

/**
 * @brief  Interface to use in class ILPInterface. enables Gurobi for solving
 *
 */
class GRBInterface
{
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

    /**
     * @brief Modifies original objective by scale and offset
     *
     * @param scale
     * @param offset
     */
    void modify_objectives(const Point &scale, const Point &offset);

    /**
     * @brief Single objective weighted sum scalarization.
     *
     * @param weighting
     * @return std::pair<SolverStatus, double>
     * First element is status of the solver, second element is value of wss
     */
    std::pair<SolverStatus, double> wss(const Point &weighting);

    /**
     * @brief Enables lexicographic solving
     *
     * @param weighting Not used in objectives, only for final wss
     * @param rel_weight attr ObjNWeight
     * @param prio attr ObjNPriority
     * @return std::pair<SolverStatus, double>
     * First element is status of the solver, second element is value of wss
     */
    std::pair<SolverStatus, double> lex_wss(const Point &weighting, const Point &rel_weight,
                                            const std::vector<int> &prio);

    /**
     * @brief Returns number of constraints
     *
     * @return int
     */
    inline int m() const
    {
        return m_m;
    }

    /**
     * @brief Returns number of variables
     *
     * @return int
     */
    inline int n() const
    {
        return m_n;
    }

    /**
     * @brief Returns number of objectives
     *
     * @return unsigned int
     */
    inline unsigned int d() const
    {
        return m_d;
    }

    /**
     * @brief Transforms unscaled point of scaled solution
     *
     * @param point Scaled solution
     * @return Point
     */
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

    /**
     * @brief For each objective: 1 if it minimized originally, -1 else
     *
     * @return const Point&
     */
    inline const Point &sense_multi() const
    {
        return m_sense_multi;
    }

    /**
     * @brief Scale factor for each objective
     *
     * @return const Point&
     */
    inline const Point &scale() const
    {
        return m_scale;
    }

    /**
     * @brief Offset for each objective
     *
     * @return const Point&
     */
    inline const Point &offset() const
    {
        return m_offset;
    }

    /**
     * @brief Value of each objective in last solution
     *
     * @param d
     * @return double
     */
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

    /**
     * @brief Name of variable j
     *
     * @param j
     * @return std::string
     */
    inline std::string var_name(int j) const
    {
        assert(j >= 0 && j < m_n);
        return m_vars_last[j].get(GRB_StringAttr_VarName);
    }

    /**
     * @brief Value of variable j in last solution
     *
     * @param j
     * @return double
     */
    inline double var_value(int j) const
    {
        assert(j >= 0 && j < m_n);
        return m_vars_last[j].get(GRB_DoubleAttr_X);
    }

    /**
     * @brief Values of all variables in last solution.
     * Returned type must allow [] operator and be self destructible
     *
     * @return std::unique_ptr<double[]>
     */
    inline std::unique_ptr<double[]> var_values()
    {
        if (m_last_is_normal)
        {
            return std::unique_ptr<double[]>{m_lp.get(GRB_DoubleAttr_X, m_vars_last.get(), m_n)};
        }
        else
        {
            return std::unique_ptr<double[]>{
                m_lex_lp.get(GRB_DoubleAttr_X, m_vars_last.get(), m_n)};
        }
    }

    /**
     * @brief Return type of var j
     *
     * @param j
     * @return VarType
     */
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

#endif
