#pragma once
/*
 * kurz_parser.h
 *
 *  Created on: 08.11.2013
 *      Author: fritz
 */

#ifndef KURZ_PARSER_H_
#define KURZ_PARSER_H_

#include <list>
#include <string>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

class KurzParser {
public:
	KurzParser(std::string filename) :
		filename_(filename) {}

	unsigned int get_graph(ogdf::Graph &graph, ogdf::NodeArray<Point *> &node_weights, ogdf::EdgeArray<Point *> &edge_weights);
	unsigned int get_graph(ogdf::Graph &graph, ogdf::EdgeArray<Point *> &edge_weights);
	unsigned int get_graph(ogdf::Graph &graph, ogdf::NodeArray<Point *> &node_weights);

private:
	std::string filename_;

};

} /* namespace mco */
#endif /* KURZ_PARSER_H_ */
