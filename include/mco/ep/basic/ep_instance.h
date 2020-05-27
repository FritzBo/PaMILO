#pragma once
/*
 * EpInstance.h
 *
 *  Created on: 12.03.2013
 *      Author: fritz
 */

#ifndef EPINSTANCE_H_
#define EPINSTANCE_H_

#include <exception>
#include <string>

#include <ogdf/basic/Graph.h>

#include <mco/basic/abstract_graph_instance.h>

namespace mco {

class Point;

class EpInstance : public AbstractGraphInstance {

public:
	EpInstance(ogdf::Graph graph,
               ogdf::EdgeArray<Point *> weights,
               unsigned int dimension,
               ogdf::node const source,
               ogdf::node const target)
    :   AbstractGraphInstance(graph, weights, dimension),
        source_(source),
        target_(target) {}

	ogdf::node source() const { return source_; }
	ogdf::node target() const { return target_; }

private:
	ogdf::node const source_;
	ogdf::node const target_;
};

}

#endif /* EPINSTANCE_H_ */
