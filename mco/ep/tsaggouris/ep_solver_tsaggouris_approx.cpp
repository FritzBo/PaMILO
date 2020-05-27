/*
 * ep_solver_tsaggouris_approx.cpp
 *
 *  Created on: 26.03.2013
 *      Author: fritz
 */

#include <mco/ep/tsaggouris/ep_solver_tsaggouris_approx.h>

#include <vector>
#include <list>
#include <cmath>
#include <stdexcept>
#include <functional>

using std::invalid_argument;
using std::vector;
using std::list;
using std::log;
using std::floor;
using std::pair;
using std::max;
using std::min;
using std::numeric_limits;
using std::function;

#include <ogdf/basic/Graph.h>

using ogdf::node;
using ogdf::edge;
using ogdf::NodeArray;
using ogdf::EdgeArray;
using ogdf::Graph;

#include <mco/ep/basic/ep_instance.h>
#include <mco/ep/martins/label.h>

namespace mco {

EpSolverTsaggourisApprox::EpSolverTsaggourisApprox(EpInstance &instance, const Point epsilon) : AbstractEpSolver(instance), epsilon_(epsilon) {
}

unsigned int EpSolverTsaggourisApprox::position(const Point &point) const {
	unsigned int pos = 0;
	for(unsigned int i = 0; i < instance().dimension() - 1; ++i)
		pos += bases_[i] * static_cast<unsigned int>(floor(log(point[i]/c_min_[i])/log(epsilon_[i])));
	return pos;
}

void EpSolverTsaggourisApprox::ExtendAndMerge(vector<const Label *> &old_Py_n, const edge e, const node n, vector<const Label *> &new_Py_n) const {
	const function<Point *(edge)> & weights = instance().weights();
	const unsigned int dimension = instance().dimension();

	for(auto label : old_Py_n) {
		if(label == nullptr)
			continue;

		const Label *new_label = new const Label(new Point(*label->point + *weights(e)), n, label);
		const Label *old_label = new_Py_n[position(*new_label->point)];
		if(old_label == nullptr || (*new_label->point)[dimension - 1] < (*old_label->point)[dimension - 1]) {
			new_Py_n[position(*new_label->point)] = new_label;
			delete old_label;
		} else
			delete new_label;
	}
}

void EpSolverTsaggourisApprox::Solve() {
	const unsigned int dimension = instance().dimension();
	const Graph &graph = instance().graph();
	const function<Point *(edge)> & weights = instance().weights();

	vector<double> c_max(dimension - 1);
	c_min_.reserve(dimension - 1);

	for(unsigned int i = 0; i < dimension - 1; ++i) {
		c_max[i] = 0;
		c_min_[i] = numeric_limits<double>::infinity();
	}

	edge e;
	forall_edges(e, graph) {
		if(e->isSelfLoop())
			continue;

		for(unsigned int i = 0; i < dimension - 1; ++i) {
			c_max[i] = max(c_max[i], (*weights(e))[i]);
			c_min_[i] = min(c_min_[i], (*weights(e))[i]);

			if(c_min_[i] == 0)
				throw invalid_argument("Edge cost 0 is not allowed.");
		}
	}

	double a;
	unsigned int size = 1;
	bases_.reserve(dimension - 1);

	for(unsigned int i = 0; i < dimension - 1; ++i) {
		a = log(graph.numberOfNodes() * (c_max[i] / c_min_[i])) / log(epsilon_[i]);
		bases_[i] = size;
		size *= static_cast<unsigned int>(floor(a));
	}

	NodeArray<vector<const Label *>> * old_Py = new NodeArray<vector<const Label *>>(graph);

	node n;
	forall_nodes(n, graph) {
		(*old_Py)[n].resize(size, nullptr);
	}

	(*old_Py)[instance().source()][0] = new Label(Point::Null(dimension), nullptr, nullptr);

	for(int i = 1; i < graph.numberOfNodes(); ++i)
		forall_nodes(n, graph) {
			auto new_Py = new NodeArray<vector<const Label *>>(*old_Py);
			forall_adj_edges(e, n) {
				if(e->target() != n)
					continue;

				ExtendAndMerge((*old_Py)[e->source()], e, e->source(), (*new_Py)[n]);
			}
			delete old_Py;
			old_Py = new NodeArray<vector<const Label *>>(*new_Py);
			delete new_Py;
		}

	list<pair<csolution_type, const Point>> target_points;
	for(auto label : (*old_Py)[instance().target()])
		if(label != nullptr)
			target_points.push_back(make_pair(list<edge>(), *label->point));

	add_solutions(target_points.begin(), target_points.end());

	delete old_Py;
}

} /* namespace mco */
