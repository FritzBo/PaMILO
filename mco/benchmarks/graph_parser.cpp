/*
 * graph_parser.cpp
 *
 *  Created on: 10.06.2013
 *      Author: fritz
 */

#include <mco/benchmarks/graph_parser.h>

#include <map>

using std::map;
using std::pair;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::node;
using ogdf::edge;
using ogdf::EdgeArray;

#include <mco/basic/point.h>

using mco::Point;

namespace mco {

void GraphParser::getGraph(Graph &graph, EdgeArray<Point *> &weights) {

	ifstream file(filename_);

	if(!file.good()) {
	    cerr << "Could not open file " << filename_ << endl;
        throw string("Could not open file ") + filename_;
    }

	int n, m;

	file >> n; //number of nodes
	file >> m; //number of variables
	int dim = 2;

	int node1_ref, node2_ref;
	node node1, node2;
	double * objective_value;
	edge e;
	map<int, node> nodes_added;

	for( int i = 0; i < m; i++ ) {
		file >> node1_ref >> node2_ref;

		if(nodes_added.count(node1_ref) == 0) {
			node1 = graph.newNode(node1_ref);
			nodes_added.insert(pair<int, node>(node1_ref, node1));
		} else
			node1 = nodes_added[node1_ref];

		if(nodes_added.count(node2_ref) == 0) {
			node2 = graph.newNode(node2_ref);
			nodes_added.insert(pair<int, node>(node2_ref, node2));
		} else
			node2 = nodes_added[node2_ref];

		objective_value = new double[dim];
		for(int j = 0; j < dim; j++) {
			file >> objective_value[j];
		}

		e = graph.newEdge(node1, node2);
		weights[e] = new Point(objective_value, dim);
		delete[] objective_value;

	}

	//close the file
	file.close();
}

}
