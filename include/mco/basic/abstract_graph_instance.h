#pragma once
/*
 * vector_weighted_graph.h
 *
 *  Created on: 25.03.2013
 *      Author: fritz
 */

#ifndef ABSTRACT_MCO_GRAPH_INSTANCE_H_
#define ABSTRACT_MCO_GRAPH_INSTANCE_H_

#include <functional>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

class AbstractGraphInstance  {
    
public:
    
	AbstractGraphInstance(const ogdf::Graph & graph,
                          const ogdf::EdgeArray<Point *> & weights,
                          unsigned int dimension)
    : graph_(graph), weights_(weights), dimension_(dimension) { }

	AbstractGraphInstance(const AbstractGraphInstance & graph) = delete;
	void operator=(const AbstractGraphInstance & graph) = delete;

	unsigned int dimension() const {
        return dimension_;
    }
    
    const ogdf::EdgeArray<Point*>& weights() const {
        return weights_;
    }
    
    const ogdf::Graph & graph() const {
        return graph_;
    }
    
private:
    
    const ogdf::Graph & graph_;
	const ogdf::EdgeArray<Point *> & weights_;
	unsigned int dimension_;

};

} /* namespace mco */

#endif /* ABSTRACT_MCO_GRAPH_INSTANCE_H_ */
