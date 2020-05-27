#pragma once
/*
 * bsssa.h
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#ifndef _EP_WEIGHTED_BS_
#define _EP_WEIGHTED_BS_

#include <mco/basic/abstract_solver.h>

namespace mco {

class EpWeightedBS : public AbstractSolver<std::list<ogdf::edge>> {
    
public:
	EpWeightedBS(double epsilon = 0)
    : epsilon_(epsilon) { }
    
	virtual void Solve(const ogdf::Graph& graph,
                       std::function<const Point*(const ogdf::edge)> costs,
                       unsigned dimension,
                       const ogdf::node source,
                       const ogdf::node target,
                       bool directed = true);
    
private:
    bool ConvexHull(const std::list<const Point *> &source1,
                    const std::list<const Point *> &source2,
                    std::list<const Point *> &nondominated_subset,
                    std::list<const Point *> &dominated_subset,
                    double epsilon);

    const double epsilon_;
    
};

}

#endif /* _EP_WEIGHTED_BS_ */