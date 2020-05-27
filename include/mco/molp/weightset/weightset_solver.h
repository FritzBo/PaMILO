#pragma once
/*
 * weightset_solver.h
 *
 *  Created on: 19.07.2013
 *      Author: fritz
 */

#ifndef WEIGHTSET_SOLVER_H_
#define WEIGHTSET_SOLVER_H_

#include <list>

#include <gurobi_c++.h>

#include <mco/molp/basic/molp_model.h>
#include <mco/basic/abstract_solver.h>

namespace mco {

class WeightSetMolpSolver : AbstractSolver<std::list<unsigned>> {
	const MolpModel &molp_model_;
	GRBEnv *grb_env_;
	GRBModel *grb_model_;
	GRBVar *x_;
	GRBConstr *c_;
	GRBConstr weight_lock_;
	std::list<GRBConstr> lock_constraints_;

	void Init();
	void set_weights(double *weights);
	void lock_weight(double *weights, double value);
	void optimize_objective(unsigned int objective, char sense);
	void lock_objective_value(unsigned int objective, double value);
	void unlock();

	bool epsilon_equal(double * point1, double * point2, unsigned int dimension, double epsilon);
	bool is_dominated(double * point1, double * point2, unsigned int dimension, double epsilon);

	void cone_combination(list<double *> vertices, double * combination, unsigned int dimension);
	double inner_product(double * vector1, double * vector2, unsigned int dimension);

	void lexopt_on_face(double *l, char sense, unsigned int variables, unsigned int objectives, double *x, double *y);

	string print_point(double * point, unsigned int dimension);

public:
	WeightSetMolpSolver() = delete;
	explicit WeightSetMolpSolver(const MolpModel &model);
	WeightSetMolpSolver(GRBEnv * grb_env, const MolpModel &model);

	void Solve();
};

}


#endif /* WEIGHTSET_SOLVER_H_ */
