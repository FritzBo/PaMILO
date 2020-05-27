/*
 * kurz_parser.cpp
 *
 *  Created on: 08.11.2013
 *      Author: fritz
 */

#include <mco/benchmarks/kurz_parser.h>

#include <vector>
#include <exception>
#include <string>
#include <istream>
#include <list>

using std::getline;
using std::exception;
using std::vector;
using std::list;

#include <ogdf/basic/Graph.h>

using ogdf::node;
using ogdf::edge;

namespace mco {

} /* namespace mco */

unsigned int mco::KurzParser::get_graph(ogdf::Graph& graph,
		ogdf::NodeArray<Point *>& node_weights,
		ogdf::EdgeArray<Point *>& edge_weights) {

	if(node_weights.graphOf() != &graph || edge_weights.graphOf() != &graph)
		return 0;

	return 0;
}

unsigned int mco::KurzParser::get_graph(ogdf::Graph& graph,
		ogdf::EdgeArray<Point *>& edge_weights) {

	if(edge_weights.graphOf() != &graph)
		return 0;

	ifstream file(filename_);

	if(!file.good())
		return 0;

	string str;
	getline(file, str);
	int number_nodes = stoi(str.substr(str.find_first_of('=', 0) + 2, str.length() - 1), 0, 10);

	vector<node> node_indices(number_nodes, nullptr);

	unsigned int dimension = 0;
	int source_id, target_id, number_attributes;
	double value;
	list<double> values;
	Point *p;
	edge e;

	for(int i = 0; i < number_nodes; ++i) {
		file >> source_id;

		if(node_indices[source_id] == nullptr)
			node_indices[source_id] = graph.newNode(source_id);

		while(file.peek() == ':') {
			file.ignore(1);
			file >> value;
			values.push_back(value);
		}
		if(!values.empty()) {
			number_attributes = values.size();
//			p = new Point(number_attributes);
			for(int j = 0; j < number_attributes; ++j) {
//				p->operator [](j) = values.front();
				values.pop_front();
			}
		}

		while(file.peek() == ' ') {
			file >> target_id;
			if(node_indices[target_id] == nullptr)
				node_indices[target_id] = graph.newNode(target_id);

			e = graph.newEdge(node_indices[source_id], node_indices[target_id]);

			values.clear();
			while(file.peek() == ':') {
				file.ignore(1);
				file >> value;
				values.push_back(value);
			}
			if(!values.empty()) {

				if(dimension == 0)
					dimension = values.size();
				else if(dimension != values.size())
					return 0;

				p = new Point(dimension);
				for(int j = 0; j < dimension; ++j) {
					p->operator [](j) = values.front();
					values.pop_front();
				}
				edge_weights[e] = p;
			}
		}
	}

	file.close();

	return dimension;
}

unsigned int mco::KurzParser::get_graph(ogdf::Graph& graph,
		ogdf::NodeArray<Point *>& node_weights) {

	if(node_weights.graphOf() != &graph)
			return 0;

	return 0;
}
