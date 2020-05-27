/*
 * prim_st_solver.cpp
 *
 *  Created on: 27.09.2013
 *      Author: fritz
 */

#include <mco/est/basic/kruskal_st_solver.h>

#include <vector>
#include <list>

using std::vector;
using std::list;

using ogdf::edge;
using ogdf::node;

namespace mco {

KruskalSTSolver::KruskalSTSolver(const ogdf::Graph &g, const ogdf::EdgeArray<double>& costs) :
	graph_(g), costs_(costs), cost_(0) {
}

KruskalSTSolver::~KruskalSTSolver() {
}

void KruskalSTSolver::Solve() {
	vector<edge> sorted_edges(graph_.numberOfEdges());

	cost_ = 0;
	spanning_tree_.clear();
	spanning_tree_.resize(graph_.numberOfNodes() - 1);

	edge e;
	forall_edges(e, graph_) {
		sorted_edges.push_back(e);
	}

	sort(sorted_edges.begin(), sorted_edges.end(), [this] (edge e1, edge e2) {return costs_[e1] < costs_[e2];});

	node n;
	forall_nodes(n, graph_) {
		make_set(n);
	}

	for(edge e: sorted_edges) {
		if(find_set(e->source()) != find_set(e->target())) {
			spanning_tree_.push_back(e);
			set_union(e->source(), e->target());
			cost_ += costs_[e];
		}
	}
}

const vector<edge> & KruskalSTSolver::spanning_tree() {
	return spanning_tree_;
}

double KruskalSTSolver::min_cost() {
	return cost_;
}

void KruskalSTSolver::make_set(ogdf::node n) {
	sets_[n] = n;
}

ogdf::node KruskalSTSolver::find_set(ogdf::node n) {
	node parent = n;
	list<node> path_nodes;

	while(sets_[parent] != parent)
		path_nodes.push_back(parent);
		parent = sets_[parent];

	for(node n: path_nodes)
		sets_[n] = parent;

	return parent;
}

void KruskalSTSolver::set_union(ogdf::node u, ogdf::node v) {
	sets_[find_set(v)] = sets_[find_set(u)];
}

} /* namespace mco */
