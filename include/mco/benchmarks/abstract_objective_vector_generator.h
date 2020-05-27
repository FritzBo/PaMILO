#pragma once
/*
 * abstract_objective_vector_generator.h
 *
 *  Created on: 12.08.2013
 *      Author: fritz
 */

#ifndef ABSTRACT_OBJECTIVE_VECTOR_GENERATOR_H_
#define ABSTRACT_OBJECTIVE_VECTOR_GENERATOR_H_

#include <list>

#include <mco/basic/point.h>

namespace mco {

class AbstractObjectiveVectorGenerator {

public:
	explicit AbstractObjectiveVectorGenerator(unsigned int dimension) : dimension_(dimension) {}

	AbstractObjectiveVectorGenerator() = delete;
	AbstractObjectiveVectorGenerator(const AbstractObjectiveVectorGenerator &) = delete;
	AbstractObjectiveVectorGenerator & operator=(AbstractObjectiveVectorGenerator &) = delete;

	virtual Point *draw_point() = 0;

	unsigned int dimension() {
		return dimension_;
	}

	virtual ~AbstractObjectiveVectorGenerator() {
		for(auto point : points_)
			delete point;
	}

protected:
	std::list<const Point *> points_;

private:
	unsigned int dimension_;
};

} /* namespace mco */
#endif /* ABSTRACT_OBJECTIVE_VECTOR_GENERATOR_H_ */
