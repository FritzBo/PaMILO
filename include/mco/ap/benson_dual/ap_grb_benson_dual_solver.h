#pragma once
/*
 * ap_benson_dual_solver.h
 *
 *  Created on: 14.10.2013
 *      Author: fritz
 */

#ifndef AP_GRB_BENSON_DUAL_SOLVER_H_
#define AP_GRB_BENSON_DUAL_SOLVER_H_

#include <memory>

#include <gurobi_c++.h>

#include <mco/ap/basic/abstract_ap_solver.h>
#include <mco/generic/benson_dual/dual_benson_scalarizer.h>

namespace mco {

template<typename OnlineVertexEnumerator>
class APGRBBensonDualSolver :
    public AbstractAPSolver {
        
public:
	APGRBBensonDualSolver(AssignmentInstance & instance, double epsilon = 1E-8);

	virtual ~APGRBBensonDualSolver();

	void Solve() {
		std::list<Point *> solutions;
		dual_benson_solver_->Calculate_solutions(solutions);

		add_solutions(solutions.begin(), solutions.end());
	}

	double get_solver_time() {
		return cycles_ / (double) CLOCKS_PER_SEC;
	}

protected:
	virtual double Solve_scalarization(Point &weights, Point &value);

private:
	GRBEnv *grb_env_;
	GRBModel *lp_model_;

	GRBVar *variables_;

	double epsilon_;

	DualBensonScalarizer<OnlineVertexEnumerator> *dual_benson_solver_;

	clock_t cycles_;
};

} /* namespace mco */
#endif /* AP_GRB_BENSON_DUAL_SOLVER_H_ */
