/*
 * pge10_parser3_ap.cpp
 *
 *  Created on: 06.11.2013
 *      Author: fritz
 */

#include <mco/benchmarks/pge10_parser3_ap.h>

#include <set>
#include <list>

using std::make_shared;
using std::set;
using std::list;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::node;
using ogdf::edge;

namespace mco {

AssignmentInstance * PGE10Parser3AP::get_instance() {
	ifstream file(filename_);

	if(!file.good())
		    cerr << "Could not open file " << filename_ << endl;

	unsigned int num_agents;

	file >> num_agents; //number of nodes
	unsigned int dim = 3;

	auto graph = new Graph();
	auto edge_cost = new EdgeArray<Point *>(*graph);
	auto agents = new set<node>();
	list<node> jobs;

	for(unsigned int i = 0; i < num_agents; ++i) {
		node agent_node = graph->newNode();
		agents->insert(agent_node);
	}

	for(unsigned int i = 0; i < num_agents; ++i) {
		node job_node = graph->newNode();
		jobs.push_back(job_node);
		for(auto agent : *agents)
			graph->newEdge(agent, job_node);
	}

	double value;
	edge e;
	for(unsigned int i = 0; i < dim; ++i)
		for(auto agent : *agents)
			forall_adj_edges(e, agent){
				if(i == 0) {
					edge_cost->operator [](e) = new Point(0.0, dim);
				}
				file >> value;
				edge_cost->operator [](e)->operator [](i) = value;
			}

	auto instance = new AssignmentInstance(*graph, *edge_cost, *agents, dim);

	//close the file
	file.close();

	return instance;
}

} /* namespace mco */
