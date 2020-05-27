/*
 * est_solver2_trees.cpp
 *
 *  Created on: 25.03.2013
 *      Author: fritz
 */

#include <mco/est/2tree/est_solver2_trees.h>

#include <list>

using std::list;

#include <ogdf/basic/Graph.h>

using ogdf::edge;
using ogdf::node;
using ogdf::EdgeArray;
using ogdf::AdjElement;

#include <mco/basic/point.h>

using mco::Point;

namespace mco {

void ESTSolver2Trees::Solve() {
/*	EdgeArray<list<Point *>> st(instance().graph());
	EdgeArray<list<Point *>> dt(instance().graph());

	list<node> degree_2_nodes;

	node v;
	forall_nodes(v, instance().graph())
		if(v->outdeg() == 2)
			degree_2_nodes.push_back(v);

	while(!degree_2_nodes.empty()) {

		v = degree_2_nodes.front();
		degree_2_nodes.pop_front();

		node neighbor[2] = {nullptr, nullptr};
		edge L[2] = {nullptr, nullptr};
		edge R[2] = {nullptr, nullptr};
		edge M[2] = {nullptr, nullptr};

		edge e;
		forall_adj_edges(e, v) {
			if(e->source() == v) {
				if(neighbor[1] == nullptr)
					neighbor[1] = e->target();
				else if(neighbor[0] == nullptr)
					neighbor[0] = e->target();

				if(e->target() == neighbor[1])
					R[0] = e;
				else
					R[1] = e;
			} else {
				if(neighbor[0] == nullptr)
					neighbor[0] = e->source();
				else if(neighbor[1] == nullptr)
					neighbor[1] = e->source();

				if(e->source() == neighbor[0])
					L[0] = e;
				else
					L[1] = e;
			}
		}

		for(int i = 0; i < 2; ++i) {
			edge M = M[i];
			edge R = R[i];
			edge L = L[i];
			node s = neighbor[1-i];
			node t = neighbor[i];
		}

	}*/
}

} /* namespace mco */
