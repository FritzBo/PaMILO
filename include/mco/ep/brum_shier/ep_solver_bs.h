#pragma once
/*
 * bsssa.h
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#ifndef BSSSA_H_
#define BSSSA_H_

#include <mco/basic/abstract_solver.h>

namespace mco {

class EpSolverBS : public AbstractSolver<std::list<ogdf::edge>> {
    
public:
	EpSolverBS(double epsilon = 0)
    : epsilon_(epsilon) { }
    
	virtual void Solve(const ogdf::Graph& graph,
                       std::function<const Point*(const ogdf::edge)> costs,
                       unsigned dimension,
                       const ogdf::node source,
                       const ogdf::node target,
                       bool directed = true);
    
private:
    const double epsilon_;
};

}

#endif /* BSSSA_H_ */
