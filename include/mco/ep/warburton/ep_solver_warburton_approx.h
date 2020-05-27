#pragma once
/*
 * ep_solver_warburton_approx.h
 *
 *  Created on: 02.04.2013
 *      Author: fritz
 */

#ifndef EP_SOLVER_WARBURTON_APPROX_H_
#define EP_SOLVER_WARBURTON_APPROX_H_

#include <mco/ep/basic/abstract_ep_solver.h>
#include <mco/ep/basic/ep_instance.h>

namespace mco {

class EpSolverWarburtonApprox : public AbstractEpSolver {

	const Point epsilon_;
	const double  theta_;
	const unsigned int processes_;

public:
	EpSolverWarburtonApprox(EpInstance &instance, const Point &epsilon, unsigned int processes = 2, double theta = 2.0);
	void Solve();

	~EpSolverWarburtonApprox() noexcept {}
};

}

#endif /* EP_SOLVER_WARBURTON_APPROX_H_ */
