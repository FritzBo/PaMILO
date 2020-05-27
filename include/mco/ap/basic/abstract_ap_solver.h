#pragma once
/*
 * AbstractAPSolver.h
 *
 *  Created on: 14.10.2013
 *      Author: fritz
 */

#ifndef ABSTRACT_AP_SOLVER_H_
#define ABSTRACT_AP_SOLVER_H_

#include <memory>

#include <mco/basic/abstract_solver.h>
#include <mco/ap/basic/ap_instance.h>

namespace mco {
    
    // FIXME

class AbstractAPSolver :
public AbstractSolver<std::list<ogdf::edge>> {
        
public:
        
	inline AbstractAPSolver(AssignmentInstance & instance)
    :   instance_(instance) {
    }
    
    virtual ~AbstractAPSolver() { }

protected:
        
	inline AssignmentInstance & instance() {
		return instance_;
	}

private:
	AssignmentInstance & instance_;

};
    
}


#endif /* ABSTRACT_AP_SOLVER_H_ */
