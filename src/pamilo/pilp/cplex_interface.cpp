/**
 * @file cplex_interface.chpp
 * @author Levin Nemesch
 * @brief
 * @date 08.02.2022
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */
#ifdef USE_CPLEX

#    include <pamilo/pilp/cplex_interface.hpp>

namespace pamilo {
CPLEXInterface::CPLEXInterface(PilpBensonArgs &args)
    : m_env()
    , m_model(m_env)
    , m_cplex(m_env)
    , m_vars(m_env)
    , m_cons(m_env)
{
    std::string output_name = args.output_name.getValue();
    if (output_name == "")
    {
        output_name = args.instance_name.getValue();
    }
    m_cplexFile = std::ofstream(output_name + "_cplex");

    m_cplex.importModel(m_model, args.instance_name.getValue().c_str(), m_og_obj, m_vars, m_cons);

    m_cplex.setParam(IloCplex::Param::MultiObjective::Display, 2);
    m_cplex.setParam(IloCplex::Param::ParamDisplay, 0);
    m_cplex.setOut(m_cplexFile);
    m_cplex.setParam(IloCplex::Param::Threads, args.solver_threads_limit.getValue() <= 1
                                                   ? 1
                                                   : args.solver_threads_limit.getValue());
    m_cplex.extract(m_model);

    m_d = m_og_obj.getNumCriteria();
    m_m = m_cplex.getNrows();
    m_n = m_cplex.getNcols();

    m_scale = Point(1.0, m_d);
    m_offset = Point(0.0, m_d);

    if (m_og_obj.getSense() == IloObjective::Sense::Maximize)
    {
        m_sense_multi = Point(-1.0, m_d);
    }
    else
    {
        m_sense_multi = Point(1.0, m_d);
        m_scaled_obj = m_og_obj;
    }

    // clean tolerances:
    IloNumExprArray objs_tmp(m_env);
    IloNumArray weights_tmp(m_env);
    IloIntArray prio_tmp(m_env);
    IloNumArray absTols_tmp(m_env);
    IloNumArray relTols_tmp(m_env);
    for (unsigned int i = 0; i < m_d; i++)
    {
        objs_tmp.add(m_og_obj.getCriterion(i) * m_sense_multi[i]);
        weights_tmp.add(1);
        prio_tmp.add(0);
        absTols_tmp.add(0);
        relTols_tmp.add(0);
    }
    m_scaled_obj = IloObjective(
        m_env, IloStaticLex(m_env, objs_tmp, weights_tmp, prio_tmp, absTols_tmp, relTols_tmp),
        IloObjective::Minimize);

    m_last_obj = m_scaled_obj;

    m_model.remove(m_og_obj);
    m_model.add(m_last_obj);
    m_cplex.extract(m_model);
}

void CPLEXInterface::modify_objectives(const Point &scale, const Point &offset)
{
    assert(scale.dimension() == m_d && offset.dimension() == m_d);

    m_scale = scale;
    m_offset = offset;

    IloNumExprArray objs_tmp(m_env);
    IloNumArray weights_tmp(m_env);
    IloIntArray prio_tmp(m_env);
    IloNumArray absTols_tmp(m_env);
    IloNumArray relTols_tmp(m_env);
    for (unsigned int i = 0; i < m_d; i++)
    {
        objs_tmp.add((m_sense_multi[i] * m_og_obj.getCriterion(i) + m_offset[i]) * m_scale[i]);
        weights_tmp.add(1);
        prio_tmp.add(0);
        absTols_tmp.add(0);
        relTols_tmp.add(0);
    }
    m_scaled_obj = IloObjective(
        m_env, IloStaticLex(m_env, objs_tmp, weights_tmp, prio_tmp, absTols_tmp, relTols_tmp),
        IloObjective::Minimize);
}

std::pair<SolverStatus, double> CPLEXInterface::wss(const Point &weighting)
{
    IloNumExpr single_obj(m_env);
    for (unsigned int i = 0; i < m_d; i++)
    {
        single_obj += weighting[i] * m_scaled_obj.getCriterion(i);
    }

    m_model.remove(m_last_obj);

    m_last_obj = IloObjective(m_env, single_obj, IloObjective::Minimize);
    m_model.add(m_last_obj);

    m_cplex.extract(m_model);

    m_cplex.solve();

    auto status = m_cplex.getStatus();

    if (status == IloAlgorithm::Status::Optimal)
    {
        return std::pair<SolverStatus, double>{SolverStatus::Success, m_cplex.getObjValue()};
    }
    else if (status == IloAlgorithm::Status::Infeasible)
    {
        return std::pair<SolverStatus, double>{SolverStatus::Infeasible, -1};
    }
    else if (status == IloAlgorithm::Status::Unbounded)
    {
        return std::pair<SolverStatus, double>{SolverStatus::Unbounded, -1};
    }
    else
    {
        return std::pair<SolverStatus, double>{SolverStatus::UnknownStatus, -1};
    }
}

std::pair<SolverStatus, double> CPLEXInterface::lex_wss(const Point &weighting,
                                                        const Point &rel_weight,
                                                        const std::vector<int> &prio)
{
    IloNumExpr single_obj(m_env);
    IloNumExprArray objs_tmp(m_env);
    IloNumArray weights_tmp(m_env);
    IloIntArray prio_tmp(m_env);
    IloNumArray absTols_tmp(m_env);
    IloNumArray relTols_tmp(m_env);

    for (unsigned int i = 0; i < m_d; i++)
    {
        objs_tmp.add(m_scaled_obj.getCriterion(i));
        weights_tmp.add(rel_weight[i]);
        prio_tmp.add(prio[i]);
        absTols_tmp.add(0);
        relTols_tmp.add(0);

        single_obj += weighting[i] * m_scaled_obj.getCriterion(i);
    }

    m_model.remove(m_last_obj);

    m_last_obj = IloObjective(
        m_env, IloStaticLex(m_env, objs_tmp, weights_tmp, prio_tmp, absTols_tmp, relTols_tmp),
        IloObjective::Minimize);
    m_model.add(m_last_obj);

    m_cplex.extract(m_model);

    m_cplex.solve();

    auto status = m_cplex.getStatus();

    if (status == IloAlgorithm::Status::Optimal)
    {
        return std::pair<SolverStatus, double>{SolverStatus::Success, m_cplex.getValue(single_obj)};
    }
    else if (status == IloAlgorithm::Status::Infeasible)
    {
        return std::pair<SolverStatus, double>{SolverStatus::Infeasible, -1};
    }
    else if (status == IloAlgorithm::Status::Unbounded)
    {
        return std::pair<SolverStatus, double>{SolverStatus::Unbounded, -1};
    }
    else
    {
        return std::pair<SolverStatus, double>{SolverStatus::UnknownStatus, -1};
    }
}

}  // namespace pamilo

#endif
