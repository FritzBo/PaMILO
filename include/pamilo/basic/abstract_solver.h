//
//  abstract_solver.h
//
//  Created on: 25.03.2013
//      Author: Fritz Bökler
//
//  This file is distributed under the terms of
//
//  the GNU General Public License v3,
//  a copy of which can be found in the file LICENCE-GPLv3.txt
//
//  OR
//
//  for academics, a MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
//

#pragma once

#include <list>
#include <iterator>

#include <pamilo/basic/point.h>

namespace pamilo {

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

