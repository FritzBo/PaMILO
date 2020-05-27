/*
 * mcap_parser.cpp
 *
 *  Created on: 11.11.2013
 *      Author: fritz
 */

#include <mco/benchmarks/mcap_parser.h>

#include <set>
#include <list>

using std::make_shared;
using std::set;
using std::list;
using std::shared_ptr;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::node;
using ogdf::edge;

namespace mco {

// FIXME
AssignmentInstance MCAPParser::
get_instance(ogdf::Graph& graph,
             ogdf::EdgeArray<Point *>& edge_costs,
             std::set<ogdf::node>& agents) {
    
	ifstream file(filename_);

	if(!file.good()) {
        cerr << "Could not open file " << filename_ << endl;
        throw string("Could not open file ") + filename_;
    }

	unsigned int num_agents, dim;

	file >> num_agents; //number of nodes
	file >> dim; 		//dimension of objective function;

	list<node> jobs;

	for(unsigned int i = 0; i < num_agents; ++i) {
		node agent_node = graph.newNode();
		agents.insert(agent_node);
	}

	for(unsigned int i = 0; i < num_agents; ++i) {
		node job_node = graph.newNode();
		jobs.push_back(job_node);
		for(auto agent : agents)
			graph.newEdge(agent, job_node);
	}

	double value;
    for(unsigned int i = 0; i < dim; ++i) {
        for(node agent : agents) {
            for(auto adj : agent->adjEdges) {
                edge e = adj->theEdge();
                
				if(i == 0) {
					edge_costs[e] = new Point(0.0, dim);
				}
				file >> value;
				edge_costs[e]->operator [](i) = value;
			}
        }
    }

	AssignmentInstance instance(graph, edge_costs, agents, dim);

	//close the file
	file.close();

	return instance;
}

} /* namespace mco */
