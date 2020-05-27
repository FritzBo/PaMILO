#pragma once
/*
 * prim_st_solver.h
 *
 *  Created on: 27.09.2013
 *      Author: fritz
 */

#ifndef PRIM_ST_SOLVER_H_
#define PRIM_ST_SOLVER_H_

#include <map>
#include <vector>

#include <ogdf/basic/Graph.h>

namespace mco {

class KruskalSTSolver {
public:
	KruskalSTSolver(const ogdf::Graph & g, const ogdf::EdgeArray<double> & costs);
	virtual ~KruskalSTSolver();

	void Solve();

	const std::vector<ogdf::edge> & spanning_tree();
	double min_cost();

private:
	const ogdf::Graph &graph_;
	const ogdf::EdgeArray<double> &costs_;

	std::vector<ogdf::edge> spanning_tree_;
	double cost_;

	ogdf::NodeArray<ogdf::node> sets_;
	void make_set(ogdf::node n);
	ogdf::node find_set(ogdf::node n);
	void set_union(ogdf::node u, ogdf::node v);
};

} /* namespace mco */
#endif /* PRIM_ST_SOLVER_H_ */
