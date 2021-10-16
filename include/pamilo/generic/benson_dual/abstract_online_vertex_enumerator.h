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

/**
 * @brief Interface for the Online Vertex Enumerator
 * 
 */
class AbstractOnlineVertexEnumerator {
public:
	AbstractOnlineVertexEnumerator() = delete;
	virtual ~AbstractOnlineVertexEnumerator() = default;

    /**
     * @brief Returns the cpu time
     * 
     * @return double 
     */
	double get_time() {
		return cycles_ / (double) CLOCKS_PER_SEC;
	}

	AbstractOnlineVertexEnumerator(unsigned int dimension, double epsilon) :
		dimension_(dimension),
		epsilon_(epsilon),
		cycles_(0) {
	}

    /**
     * @brief Returns wether any not enumerated vertex exists
     * 
     * @return bool
     */
    virtual inline bool has_next() = 0;

    /**
     * @brief Returns a previously not enumerated vertex
     * 
     * @return Point* 
     */
    virtual inline Point * next_vertex() = 0;

    /**
     * @brief Adds a hyperplane
     * 
     * @param vertex A point the hyperplane is anchored in
     * @param normal The normal of the hyperplane
     * @param rhs 
     */
    virtual void add_hyperplane(Point &vertex, Point &normal, double rhs) = 0;

    /**
     * @brief Returns number of hyperplanes
     * 
     * @return unsigned int 
     */
    virtual unsigned int number_of_hyperplanes() = 0;


	virtual double getDistance(Point &vertex, Point &normal, double rhs) = 0;

protected:

    /**
     * @brief dimension of the space
     * 
     */
	const unsigned int dimension_;

    /**
     * @brief epsilon for floating point comparisons
     * 
     */
	const double epsilon_;

    /**
     * @briefCPU time spent on vertex enumeration
     * 
     */
	clock_t cycles_;
};
}

