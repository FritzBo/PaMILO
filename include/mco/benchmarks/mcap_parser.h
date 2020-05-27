#pragma once
/*
 * mcap_parser.h
 *
 *  Created on: 11.11.2013
 *      Author: fritz
 */

#ifndef MCAP_PARSER_H_
#define MCAP_PARSER_H_

#include <string>
#include <memory>

#include <mco/ap/basic/ap_instance.h>

namespace mco {

class MCAPParser {
public:
	MCAPParser() = delete;

	MCAPParser(std::string filename)
    :   filename_(filename) {}

	AssignmentInstance get_instance(ogdf::Graph& graph,
                                    ogdf::EdgeArray<Point *>& edge_array,
                                    std::set<ogdf::node>& agents);

private:
	std::string filename_;
};

} /* namespace mco */
#endif /* MCAP_PARSER_H_ */
