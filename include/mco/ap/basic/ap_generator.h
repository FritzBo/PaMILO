#pragma once
/*
 * assignment_generator.h
 *
 *  Created on: 19.08.2013
 *      Author: fritz
 */

#ifndef ASSIGNMENT_GENERATOR_H_
#define ASSIGNMENT_GENERATOR_H_

#include <mco/benchmarks/abstract_objective_vector_generator.h>

#include <list>

#include <mco/ap/basic/ap_instance.h>

namespace mco {

class AssignmentGenerator {
public:
	AssignmentGenerator(unsigned int num_agents,
                        double edge_probability,
                        AbstractObjectiveVectorGenerator &objective_generator)
    
    :   objective_generator_(objective_generator),
        edge_probability_(edge_probability),
        num_agents_(num_agents) {
            
    }

	AssignmentInstance * create_instance();

private:

	AbstractObjectiveVectorGenerator &objective_generator_;
	double edge_probability_;
	unsigned int num_agents_;

};

} /* namespace mco */
#endif /* ASSIGNMENT_GENERATOR_H_ */
