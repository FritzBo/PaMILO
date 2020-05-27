/*
 * assignment_generator.cpp
 *
 *  Created on: 19.08.2013
 *      Author: fritz
 */

#include <mco/ap/basic/ap_generator.h>

#include <list>
#include <memory>
#include <set>

using std::list;
using std::set;
using std::make_shared;
using std::shared_ptr;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::node;
using ogdf::edge;
using ogdf::EdgeArray;
using ogdf::randomDouble;

namespace mco {

AssignmentInstance * AssignmentGenerator::create_instance() {

	auto graph = new Graph();
	auto weights = new EdgeArray<Point *>(*graph);
	list<node> agents;
	list<node> non_agents;

	for(unsigned int i = 0; i < num_agents_; ++i)
		agents.push_back(graph->newNode());

	for(unsigned int i = 0; i < num_agents_; ++i)
		non_agents.push_back(graph->newNode());

	for(node agent : agents)
		for(node non_agent : non_agents)
			if(randomDouble(0,1) < edge_probability_) {
				edge e1 = graph->newEdge(agent, non_agent);

				(*weights)[e1] = objective_generator_.draw_point();
			}


	auto agent_set = new set<node>(agents.cbegin(), agents.cend());
	auto assignment_instance = new AssignmentInstance(*graph, *weights, *agent_set, objective_generator_.dimension());

	return assignment_instance;

}

} /* namespace mco */
