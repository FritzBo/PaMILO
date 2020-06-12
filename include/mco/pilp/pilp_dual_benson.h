//
//  pilp_dual_benson.h
//  mco
//
//  Created by Fritz Bökler and Mirko H. Wagner on 28.05.20.
//
//

#ifndef __mco__pilp_dual_benson__
#define __mco__pilp_dual_benson__

#include <functional>

#include <mco/pilp/coin.h>

#include <mco/basic/point.h>
#include <mco/basic/abstract_solver.h>
#include <mco/basic/weight_function_adaptors.h>
#include <mco/generic/benson_dual/dual_benson_scalarizer.h>
#include <mco/generic/benson_dual/ove_fp_v2.h>

namespace mco {

class ILPSolverAdaptor {
public:
    ILPSolverAdaptor(ILP &ilp) : ilp_(ilp) {}
    
    inline double operator()(const Point& weighting,
                             Point& value);
    
    inline ~ILPSolverAdaptor();

private:
	ILP &ilp_;

    std::set<Point*, LexPointComparator> known_points_;
    
};
    
template<typename OnlineVertexEnumerator = GraphlessOVE>
class PilpDualBensonSolver : public AbstractSolver<std::list<ogdf::edge>> {
public:
    PilpDualBensonSolver(double epsilon = 1E-8)
    :   epsilon_(epsilon) {}
    
    void Solve(ILP &ilp);
    
private:
    double epsilon_;
    
};
    
// in: weigthing
// out: value (of opt sol for weighting)
// given (from constructor): ilp
// return: weighted sum of value
inline double ILPSolverAdaptor::
operator()(const Point& weighting,
           Point& value) {

	// solve ilp with weighting
	// TODO: not sure if this initialization works
	double weightedObj[ilp_.osi.getNumCols()] = {0};
	for(int i = 0; i < ilp_.dimension; i++) {
		for(auto const& [index, coeff] : ilp_.obj[i]) {
			weightedObj[index] += coeff;
		}
	}
	ilp_.osi.setObjective(weightedObj);

	// solve
	ilp_.osi.branchAndBound();


	// ilp stays as is, but obj is done by weighting
	int n = ilp_.osi.getNumCols();
	const double *sol;
	sol = ilp_.osi.getColSolution();
    
	for(int i = 0; i < ilp_.dimension; i++) {
		value[i] = 0;
		// TODO: go over obj not sol
		for(int j = 0; j < n; j++) {
			if(ilp_.obj[i].find(j) != ilp_.obj[i].end()) {
				value[i] += sol[j] * ilp_.obj[i][j];
			}
		}
	}

    if(known_points_.count(&value) == 0) {
        known_points_.insert(new Point(value));
    }
    
    return ilp_.osi.getObjValue();
}
    
ILPSolverAdaptor::~ILPSolverAdaptor() {
    for(auto point : known_points_) {
        delete point;
    }
}
 
template<typename OnlineVertexEnumerator>
inline void PilpDualBensonSolver<OnlineVertexEnumerator>::
Solve(ILP &ilp) {
    
    std::list<Point *> frontier;
    
    DualBensonScalarizer<OnlineVertexEnumerator>
    dual_benson_solver(ILPSolverAdaptor(ilp),
                       ilp.dimension,
                       epsilon_);
    
    dual_benson_solver.Calculate_solutions(frontier);
    
    std::list<std::pair<std::list<ogdf::edge>, Point>> solutions;
    
    for(auto point : frontier) {
        solutions.push_back(make_pair(std::list<ogdf::edge>(), *point));
    }
                            
    add_solutions(solutions.begin(), solutions.end());
    
}
    
}

#endif /* defined(__mco__pilp_dual_benson__) */
