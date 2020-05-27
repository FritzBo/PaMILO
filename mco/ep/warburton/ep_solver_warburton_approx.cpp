/*
 * warburton.cpp
 *
 *  Created on: 02.04.2013
 *      Author: fritz
 */

#include <vector>
#include <limits>
#include <cmath>
#include <functional>

using std::vector;
using std::min;
using std::max;
using std::log;
using std::pow;
using std::ceil;
using std::floor;
using std::numeric_limits;
using std::function;

#include <ogdf/basic/Graph.h>
#include <ogdf/graphalg/ShortestPathWithBFM.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;
using ogdf::edge;
using ogdf::ShortestPathWithBFM;

#include <mco/ep/warburton/ep_solver_warburton_approx.h>

namespace mco {

EpSolverWarburtonApprox::EpSolverWarburtonApprox(EpInstance &instance, const Point &epsilon, unsigned int processes, double theta) : AbstractEpSolver(instance), epsilon_(Point(epsilon)), theta_(theta), processes_(processes) {

}

void EpSolverWarburtonApprox::Solve() {
	const unsigned int dimension = instance().dimension();
	const unsigned int number_nodes = instance().graph().numberOfNodes();
	const Graph &graph = instance().graph();
	const function<Point *(edge)> & weights = instance().weights();
	const node source = instance().source();
	const node target = instance().target();

	vector<EdgeArray<int>> weight_functions(dimension - 1, EdgeArray<int>(graph, 0));
	NodeArray<int> distances(graph);
	NodeArray<edge> predecessor(graph);
	vector<double> min_e(dimension - 1, numeric_limits<double>::infinity());
	vector<double> max_e(dimension - 1, - numeric_limits<double>::infinity());
	vector<int> ub(dimension - 1);
	vector<int> lb(dimension - 1);
	vector<int> label_limits(dimension);
	edge e;
	double weight, d;
	int max_i, min_i;
	ShortestPathWithBFM shortest_path_module;

	// Computing the bounds
	for(unsigned int k = 0; k < dimension - 1; ++k) {

		forall_edges(e, graph) {
			weight = (*weights(e))[k];
			min_e[k] = min(min_e[k], weight);
			max_e[k] = max(max_e[k], weight);
		}

		label_limits[k] = static_cast<int>(ceil((number_nodes - 1) * theta_ / epsilon_[k]));

		max_i = static_cast<unsigned int>(ceil(log(min(max_e[k] * (number_nodes - 1), max_e[k] * (number_nodes - 1) / epsilon_[k] ) ) / log(theta_) ) + 1);
		min_i = static_cast<unsigned int>(log(max(1.0, min_e[k] * (number_nodes - 1) / epsilon_[k] ) ) / log(theta_) );

		for(int i = min_i; i < max_i; ++i) {
			d = static_cast<double>(pow(2, max_i - i));
			forall_edges(e, graph)
				weight_functions[k][e] = static_cast<int>(floor((*weights(e))[k] * (number_nodes - 1) / (epsilon_[k] * d) ));

			// TODO: Error Message / Exception
			if(!shortest_path_module.call(graph, source, weight_functions[k], distances, predecessor))
				cout << "Kein Weg gefunden!" << endl;

			if(distances[target] > (number_nodes - 1) * theta_ / epsilon_[k] || max_i - i == 1) {
				lb[k] = max_i - i + 1;
				ub[k] = static_cast<int>(ceil(log(min(max_e[k] * (number_nodes - 1), max_e[k] * (number_nodes - 1) / epsilon_[k]) / log(theta_)) + 1));
				break;
			}
		}
	}

	label_limits[dimension - 1] = numeric_limits<int>::infinity();

	cout << "bounds:" << endl;
	for(unsigned int k = 0; k < dimension - 1; ++k)
		cout << "k: " << k << ", lb: " << lb[k] << ", ub: " << ub[k] << ", limit: " << label_limits[k] << endl;

}

} // namespace mco
