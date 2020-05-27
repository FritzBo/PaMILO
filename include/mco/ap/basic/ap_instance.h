#pragma once
/*
 * assignment_instance.h
 *
 *  Created on: 19.08.2013
 *      Author: fritz
 */

#ifndef AP_INSTANCE_H_
#define AP_INSTANCE_H_

#include <set>
#include <memory>

#include <ogdf/basic/Graph.h>

#include <mco/basic/abstract_graph_instance.h>

namespace mco {

class AssignmentInstance : public AbstractGraphInstance {

public:
	AssignmentInstance(ogdf::Graph & graph,
                       ogdf::EdgeArray<Point *> & weights,
                       std::set<ogdf::node> & agents,
                       unsigned dimension)
    :   AbstractGraphInstance(graph, weights, dimension),
        agents_(agents) {
    }
    
    AssignmentInstance(const AssignmentInstance& instance)
    :   AbstractGraphInstance(instance.graph(), instance.weights(), instance.dimension()), agents_(instance.agents()) {
    }

	const std::set<ogdf::node>& agents() const {
		return agents_;
	}

private:
	const std::set<ogdf::node> & agents_;
};

} /* namespace mco */
#endif /* AP_INSTANCE_H_ */
