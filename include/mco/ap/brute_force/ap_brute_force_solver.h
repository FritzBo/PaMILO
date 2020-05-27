#pragma once
/*
 * ap_brute_force_solver.h
 *
 *  Created on: 06.11.2013
 *      Author: fritz
 */

#ifndef AP_BRUTE_FORCE_SOLVER_H_
#define AP_BRUTE_FORCE_SOLVER_H_

#include <vector>
#include <set>

#include <ogdf/basic/Graph.h>

#include <mco/ap/basic/abstract_ap_solver.h>

namespace mco {

class APBruteForceSolver : public AbstractAPSolver {
public:
	APBruteForceSolver(AssignmentInstance & instance) :
		AbstractAPSolver(instance) {}

	void Solve();

private:
	void recursive_find(unsigned int agent_index,
			std::set<ogdf::node> &jobs,
			std::list<std::list<ogdf::edge>> &matchings,
			std::list<Point> &costs,
			std::list<ogdf::edge> &current_matching,
			Point &current_cost);

	std::vector<ogdf::node> agent_list;
};

} /* namespace mco */
#endif /* AP_BRUTE_FORCE_SOLVER_H_ */
