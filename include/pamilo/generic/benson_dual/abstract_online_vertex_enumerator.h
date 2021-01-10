//
//  abstract_vertex_enumerator.h
//
//  Created on: 08.12.2013
//      Author: Fritz Bökler
//
//  This file is distributed for academics only
//  under the terms of an MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
//

#pragma once

#include <vector>
#include <map>
#include <queue>
#include <list>
#include <ctime>

#include <pamilo/basic/point.h>

namespace pamilo {

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

    virtual inline bool has_next() = 0;
	virtual inline void print_ex_points() = 0;

    virtual inline Point * next_vertex() = 0;

    virtual void add_hyperplane(Point &vertex, Point &normal, double rhs) = 0;

    virtual unsigned int number_of_hyperplanes() = 0;

	virtual double getDistance(Point &vertex, Point &normal, double rhs) = 0;

protected:

	const unsigned int dimension_;
	const double epsilon_;

	clock_t cycles_;
};
}

