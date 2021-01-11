//
//  online_vertex_enumerator.h
//
//  Created on: 30.09.2013
//      Author: fritz
//
//  This file is distributed for academics only
//  under the terms of an MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
//

#pragma once

#include <list>
#include <iostream>

#include <setoper.h>
#include <cdd.h>

#include <pamilo/basic/point.h>
#include <pamilo/generic/benson_dual/abstract_online_vertex_enumerator.h>

namespace pamilo {

class OnlineVertexEnumeratorCDD : AbstractOnlineVertexEnumerator {
public:
	OnlineVertexEnumeratorCDD() = delete;
	virtual ~OnlineVertexEnumeratorCDD();

	bool has_next() override;
	Point * next_vertex() override;
	void add_hyperplane(Point& vertex, Point& normal, double rhs) override;
	double getDistance(Point &vertex, Point &normal, double rhs) override;

	unsigned int number_of_hyperplanes() {
		return number_hyperplanes_;
	}

	double get_time() {
#ifdef NDEBUG
		std::cout << "End cycles: " << cycles_ << std::endl;
		std::cout << "Time: " << (cycles_/(double) CLOCKS_PER_SEC) << std::endl;
#endif
		return cycles_/(double) CLOCKS_PER_SEC;
	}

	OnlineVertexEnumeratorCDD(Point &initial_value, unsigned int dimension, double epsilon);

protected:
	int number_hyperplanes_;

	std::list<Point> unprocessed_vertices_;

	dd_MatrixPtr h_representation_;
};
}
