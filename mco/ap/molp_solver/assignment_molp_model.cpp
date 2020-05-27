/*
 * assignment_molp_solver.cpp
 *
 *  Created on: 19.08.2013
 *      Author: fritz
 */

#include <set>
#include <memory>
#include <functional>

using std::set;
using std::shared_ptr;
using std::function;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::edge;
using ogdf::node;
using ogdf::EdgeArray;

#include <mco/basic/point.h>
#include <mco/ap/molp_solver/assignment_molp_model.h>

namespace mco {

AssignmentMolpModel::AssignmentMolpModel(shared_ptr<AssignmentInstance> instance) :
	instance_(instance) {

	const Graph &graph = instance->graph();
	const set<node> &agents = instance->agents();
    const function<Point *(edge)> & weights = instance->weights();

	int number_of_agents = graph.numberOfNodes() / 2;
	int number_of_edges = graph.numberOfEdges() / 2;

	int number_of_constraints = number_of_agents * 4 + number_of_edges;
	int number_of_variables = number_of_edges;
	unsigned int number_of_objectives = weights(graph.firstEdge())->dimension();

	int size_of_A = number_of_constraints * number_of_variables;
	int size_of_C = number_of_variables * number_of_objectives;
	int size_of_b = number_of_constraints;

	int start_index_agent_negative = number_of_variables * number_of_agents;
	int start_index_work_positive = 2 * number_of_variables * number_of_agents;
	int start_index_work_negative = 3 * number_of_variables * number_of_agents;

	A_ = new double[size_of_A];
	C_ = new double[size_of_C];
	b_ = new double[size_of_b];

	for(int i = 0; i < size_of_A; ++i)
		A_[i] = 0;


	edge e;
	forall_edges(e, graph) {
		if(agents.count(e->source()) > 0) {
			A_[e->index() / 2 + number_of_variables * e->source()->index()] = 1;
			A_[start_index_agent_negative + e->index() / 2 + number_of_variables * e->source()->index()] = -1;

			A_[start_index_work_positive + e->index() / 2 + number_of_variables * (e->target()->index() - number_of_agents)] = 1;
			A_[start_index_work_negative + e->index() / 2 + number_of_variables * (e->target()->index() - number_of_agents)] = -1;
		} else {
			A_[e->index() / 2 + number_of_variables * e->target()->index()] = 1;
			A_[start_index_agent_negative + e->index() / 2 + number_of_variables * e->target()->index()] = -1;

			A_[start_index_work_positive + e->index() / 2 + number_of_variables * (e->source()->index() - number_of_agents)] = 1;
			A_[start_index_work_negative + e->index() / 2 + number_of_variables * (e->source()->index() - number_of_agents)] = -1;
		}

		const Point * weight = weights(e);
		for(unsigned int i = 0; i < number_of_objectives; ++i)
			C_[number_of_variables * i + e->index() / 2] = (*weight)[i];
	}

	for(int i = 0; i < number_of_variables; ++i)
		A_[4 * number_of_agents * number_of_variables + i + number_of_variables * i] = 1;

	for(int i = 0; i < number_of_agents; ++i)
		b_[i] = 1;

	for(int i = number_of_agents; i < 2 * number_of_agents; ++i)
		b_[i] = -1;

	for(int i = 2 * number_of_agents; i < 3 * number_of_agents; ++i)
		b_[i] = 1;

	for(int i = 3 * number_of_agents; i < 4 * number_of_agents; ++i)
		b_[i] = -1;

	for(int i = 4 * number_of_agents; i < size_of_b; ++i)
		b_[i] = 0;

	objectives_ = number_of_objectives;
	variables_ = number_of_variables;
	constraints_ = number_of_constraints;

}

} /* namespace mco */
