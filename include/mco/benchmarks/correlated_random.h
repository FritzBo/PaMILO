#pragma once
/*
 * correlated_random.h
 *
 *  Created on: 12.04.2013
 *      Author: fritz
 */

#ifndef CORRELATED_RANDOM_H_
#define CORRELATED_RANDOM_H_

#include <mco/benchmarks/abstract_objective_vector_generator.h>

namespace mco {

class CorrelatedObjectiveGenerator
    : public AbstractObjectiveVectorGenerator {
        
        
	double lower_bound_;
	double upper_bound_;
	double rho_;

public:
	CorrelatedObjectiveGenerator(double lower_bound,
                                 double upper_bound,
                                 double rho)
    
    : AbstractObjectiveVectorGenerator(2),
    lower_bound_(lower_bound),
    upper_bound_(upper_bound),
    rho_(rho) {}

	virtual Point * draw_point();
};

} /* namespace mco */

#endif /* CORRELATED_RANDOM_H_ */
