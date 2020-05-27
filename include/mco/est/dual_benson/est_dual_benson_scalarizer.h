#pragma once
/*
 * est_dual_benson_scalarizer.h
 *
 *  Created on: 27.09.2013
 *      Author: fritz
 */

#ifndef EST_DUAL_BENSON_SCALARIZER_H_
#define EST_DUAL_BENSON_SCALARIZER_H_

#include <vector>
#include <list>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>
#include <mco/basic/abstract_graph_instance.h>
#include <mco/est/basic/abstract_est_solver.h>
#include <mco/generic/benson_dual/dual_benson_scalarizer.h>

namespace mco {

template<class OnlineVertexEnumerator>
class ESTDualBensonScalarizer : public AbstractESTSolver {
public:
	ESTDualBensonScalarizer(AbstractGraphInstance &graph, double epsilon = 1E-6);
	virtual ~ESTDualBensonScalarizer();

	virtual void Solve() {
		std::list<Point *> solutions;
		benson_scalarizer_->Calculate_solutions(solutions);

		add_solutions(solutions.begin(), solutions.end());
	}

	double solver_time() {
		return cycles_ / (double) CLOCKS_PER_SEC;
	}

	double vertex_enumeration_time() {
		return benson_scalarizer_->vertex_enumeration_time();
	}

	int number_vertices() {
		return benson_scalarizer_->number_vertices();
	}

	int number_facets() {
		return benson_scalarizer_->number_facets();
	}


private:
	double epsilon_;
	DualBensonScalarizer<OnlineVertexEnumerator> *benson_scalarizer_;

	virtual double Solve_scalarization(Point &weights, Point &value);

	ogdf::NodeArray<ogdf::node> parents_;

	double kruskal_solver(std::vector<ogdf::edge> &sorted_edges, ogdf::EdgeArray<double> costs, Point &value);
	void make_set(ogdf::node);
	ogdf::node find_set(ogdf::node);
	void set_union(ogdf::node, ogdf::node);

	clock_t cycles_;

};

} /* namespace mco */
#endif /* EST_DUAL_BENSON_SCALARIZER_H_ */
