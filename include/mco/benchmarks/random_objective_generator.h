#pragma once
/*
 * random_objective_generator.h
 *
 *  Created on: 12.08.2013
 *      Author: fritz
 */

#ifndef RANDOM_OBJECTIVE_GENERATOR_H_
#define RANDOM_OBJECTIVE_GENERATOR_H_

#include "abstract_objective_vector_generator.h"

namespace mco {

class RandomObjectiveGenerator: public mco::AbstractObjectiveVectorGenerator {
	double lower_bound_;
	double upper_bound_;

public:
	RandomObjectiveGenerator(unsigned int dimension, double lower_bound, double upper_bound) : AbstractObjectiveVectorGenerator(dimension), lower_bound_(lower_bound), upper_bound_(upper_bound) {}

	virtual Point * draw_point();

	virtual ~RandomObjectiveGenerator();
};

} /* namespace mco */
#endif /* RANDOM_OBJECTIVE_GENERATOR_H_ */
