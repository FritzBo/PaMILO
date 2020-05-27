#pragma once
/*
 * abstract_vertex_enumerator.h
 *
 *  Created on: 08.12.2013
 *      Author: fritz
 */


#ifndef ABSTRACT_ONLINE_VERTEX_ENUMERATOR_H_
#define ABSTRACT_ONLINE_VERTEX_ENUMERATOR_H_


#include <vector>
#include <map>
#include <queue>
#include <list>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

class AbstractOnlineVertexEnumerator {
public:
	AbstractOnlineVertexEnumerator() = delete;
	virtual ~AbstractOnlineVertexEnumerator() = default;

	double get_time() {
		return cycles_ / (double) CLOCKS_PER_SEC;
	}

	AbstractOnlineVertexEnumerator(unsigned int dimension, double epsilon) :
		dimension_(dimension),
		epsilon_(epsilon),
		cycles_(0) {
	}

protected:

	const unsigned int dimension_;
	const double epsilon_;

	clock_t cycles_;
};

} /* namespace mco */

#endif /* ABSTRACT_VERTEX_ENUMERATOR_H_ */
