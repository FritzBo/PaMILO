#pragma once
/*
 * ep_solver_tsaggouris_approx.h
 *
 *  Created on: 26.03.2013
 *      Author: fritz
 */

#ifndef EP_SOLVER_TSAGGOURIS_APPROX_H_
#define EP_SOLVER_TSAGGOURIS_APPROX_H_

#include <vector>

#include <mco/basic/point.h>
#include <mco/ep/basic/abstract_ep_solver.h>
#include <mco/ep/martins/label.h>

namespace mco {

class EpSolverTsaggourisApprox: public mco::AbstractEpSolver {
	const Point epsilon_;

	unsigned int position(const Point &) const;
	void ExtendAndMerge(std::vector<const Label *> &, const ogdf::edge, const ogdf::node, std::vector<const Label *> &) const;

	std::vector<double> c_min_;
	std::vector<double> bases_;

public:
	EpSolverTsaggourisApprox(EpInstance &instance, const Point epsilon);
	void Solve();

	~EpSolverTsaggourisApprox() noexcept {}
};

} /* namespace mco */
#endif /* EP_SOLVER_TSAGGOURIS_APPROX_H_ */
