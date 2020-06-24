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
//#include <mco/generic/benson_dual/ove_fp_v2.h>
#include <mco/generic/benson_dual/ove_cdd.h>

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
    
template<typename OnlineVertexEnumerator = OnlineVertexEnumeratorCDD>
class PilpDualBensonSolver : public AbstractSolver<std::list<std::string>> {
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
	int n = ilp_.osi.getNumCols();

	// solve ilp with weighting
	//std::cout << "\nweighting: ";
	for(int i = 0; i < ilp_.dimension; i++) {
		ilp_.osi.setObjCoeff(ilp_.obj[i], weighting[ilp_.obj[i]]);
		//std::cout << weighting[i] << ", ";
	}
	//std::cout << std::endl;
	const double *coefs = ilp_.osi.getObjCoefficients();

	// solve
	ilp_.osi.branchAndBound();


	// ilp stays as is, but obj is done by weighting
	const double *sol = ilp_.osi.getColSolution();

	//std::cout << "sol:\nobj:(" << sol[0];
	//for(int i = 1; i < ilp_.dimension; i++) {
	//	std::cout << "," << sol[i];
	//}
	//std::cout << ")\nactivated vars:";
	//for(int i = ilp_.dimension; i < n; i++) {
	//	if(sol[i] == 1) {
	//		//std::cout << "(" << (i-3)%40 << "," << (i-3)/40 << ")-" << i << ", ";
	//		std::cout << i << ", ";
	//	}
	//}
	//std::cout << "\n\n";
    
	for(int i = 0; i < ilp_.dimension; i++) {
		value[i] = sol[ilp_.obj[i]];
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
    
    std::list<std::pair<std::list<std::string>, Point>> solutions;
    
    for(auto point : frontier) {
        solutions.push_back(make_pair(std::list<std::string>(), *point));
    }
                            
    add_solutions(solutions.begin(), solutions.end());
    
	std::cout << "ve time: " << dual_benson_solver.vertex_enumeration_time() << "\n";
}
    
}

#endif /* defined(__mco__pilp_dual_benson__) */
