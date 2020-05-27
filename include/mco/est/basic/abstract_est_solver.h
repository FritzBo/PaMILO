#pragma once
/*
 * abstract_est_solver.h
 *
 *  Created on: 25.03.2013
 *      Author: fritz
 */

#ifndef ABSTRACT_EST_SOLVER_H_
#define ABSTRACT_EST_SOLVER_H_

#include <list>

#include <mco/basic/abstract_solver.h>
#include <mco/basic/abstract_graph_instance.h>

namespace mco {

class AbstractESTSolver : public AbstractSolver<std::list<ogdf::edge>> {
	AbstractGraphInstance & instance_;
	std::list<Point *> solutions_;

protected:

	AbstractGraphInstance & instance() {
		return instance_;
	}

public:
	AbstractESTSolver() = delete;
	AbstractESTSolver(AbstractGraphInstance & instance) : instance_(instance) {}

};

} /* namespace mco */
#endif /* ABSTRACT_EST_SOLVER_H_ */
