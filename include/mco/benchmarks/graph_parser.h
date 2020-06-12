#pragma once
/*
 * graph_parser.h
 *
 *  Created on: 10.06.2013
 *      Author: fritz
 */

#ifndef GRAPH_PARSER_H_
#define GRAPH_PARSER_H_

#include <string>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

class GraphParser {

	std::string filename_;

public:
	GraphParser(std::string filename) :
		filename_(filename) { }

	void getGraph(ogdf::Graph &graph, ogdf::EdgeArray<mco::Point *> &weights);
};

}

#endif /* GRAPH_PARSER_H_ */
