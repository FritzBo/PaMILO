/*
 * graph_generators.cpp
 *
 *  Created on: 28.04.2013
 *      Author: fritz
 */

#include <mco/benchmarks/graph_generators.h>

#include <vector>
using std::vector;

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>
using ogdf::Graph;
using ogdf::edge;
using ogdf::node;
using ogdf::randomNumber;
using ogdf::completeGraph;

namespace mco {


// TODO: Finish
void KTreeGenerator(Graph & graph, unsigned int number_of_nodes, unsigned int k) {
	if(number_of_nodes == 0)
		return;

	if(number_of_nodes <= k) {
		completeGraph(graph, number_of_nodes);
		return;
	}

	completeGraph(graph, k);

	vector<vector<node>> candidate_cliques;

	node n;
	vector<node> initial_clique(k);
	forall_nodes(n, graph)
		initial_clique.push_back(n);
	candidate_cliques.push_back(initial_clique);

	vector<vector<node>>::iterator iter;
	vector<node> *c;
	node v;
	for(unsigned int i = 0; i < number_of_nodes - 2; ++i) {
		iter = candidate_cliques.begin();

		for(int j = 0; j < randomNumber(0, candidate_cliques.size() - 1); ++j)
			iter++;

		c = &*iter;

		candidate_cliques.erase(iter);

		v = graph.newNode();

		for(auto node: *c) {
			graph.newEdge(v, node);
			graph.newEdge(node, v);
		}
	}
}

}
