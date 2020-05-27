#pragma once
/*
 * abstract_solver.h
 *
 *  Created on: 25.03.2013
 *      Author: fritz
 */

#ifndef ABSTRACT_SOLVER_H_
#define ABSTRACT_SOLVER_H_

#include <list>
#include <iterator>

#include <mco/basic/point.h>
#include <mco/basic/abstract_graph_instance.h>

namespace mco {

template<class T>
class AbstractSolver {

public:
    
    using solution_type = T;
    using csolution_type = const T;
    
	AbstractSolver() = default;
	AbstractSolver(const AbstractSolver &) = delete;

	AbstractSolver & operator=(const AbstractSolver&) = delete;

	const std::list<std::pair<csolution_type, const Point>> & solutions() const {
        return solutions_;
    }
    
    virtual ~AbstractSolver() = default;

protected:

	void add_solution(csolution_type solution, const Point value) {
		solutions_.push_back(std::make_pair(solution, value));
	}

	template<class InputIterator>
	void add_solutions(InputIterator begin, InputIterator end) {
		solutions_.insert(solutions_.end(), begin, end);
	}

	void reset_solutions() {
		solutions_.clear();
	}

private:
    
	std::list<std::pair<csolution_type, const Point>> solutions_;

};
    
}

#endif /* ABSTRACT_SOLVER_H_ */
