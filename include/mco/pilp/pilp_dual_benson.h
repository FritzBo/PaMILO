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
#include <mco/generic/benson_dual/dual_benson_scalarizer.h>

//#define CDD_OVE
#ifdef CDD_OVE
#include <mco/generic/benson_dual/ove_cdd.h>
#else
#include <mco/generic/benson_dual/ove_fp_v2.h>
#endif

namespace mco {

class ILPSolverAdaptor {
public:
    ILPSolverAdaptor(ILP &ilp, double& time) : ilp_(ilp), time_(time) {}
    
    inline double operator()(const Point& weighting,
                             Point& value);
    
    inline ~ILPSolverAdaptor();

private:
	ILP &ilp_;

    std::set<Point*, LexPointComparator> known_points_;

	double &time_;
    
};
    
#ifdef CDD_OVE
template<typename OnlineVertexEnumerator = OnlineVertexEnumeratorCDD>
#else
template<typename OnlineVertexEnumerator = GraphlessOVE>
#endif
class PilpDualBensonSolver : public AbstractSolver<std::list<std::string>> {
public:
    PilpDualBensonSolver(double epsilon = 1E-6)
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
	for(int i = 0; i < ilp_.dimension; i++) {
		ilp_.osi.setObjCoeff(ilp_.obj[i], weighting[ilp_.obj[i]]);
	}
	const double *coefs = ilp_.osi.getObjCoefficients();

	// solve
	clock_t start = clock();
	ilp_.osi.writeLp("asdf.lp");
	ilp_.osi.branchAndBound();
	time_ += (clock() - start) / double(CLOCKS_PER_SEC);

	// ilp stays as is, but obj is done by weighting
	const double *sol = ilp_.osi.getColSolution();

	// test for integrality of solution
//	for(int i = 0; i < n; i++) {
//		if(sol[i] != int(sol[i])) {
//			std::cout << ilp_.osi.getColName(i) << " " << sol[i] << " alarm\n";
//		}
//	}

	for(int i = 0; i < ilp_.dimension; i++) {
		value[i] = sol[ilp_.obj[i]];
	}

    if(known_points_.count(&value) == 0) {
        known_points_.insert(new Point(value));
    }

	std::cout << ilp_.osi.getObjValue() << " obj\n";
    
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
	double solver_time = 0;
    
    DualBensonScalarizer<OnlineVertexEnumerator>
    dual_benson_solver(ILPSolverAdaptor(ilp,solver_time),
                       ilp.dimension,
                       epsilon_);
    
    dual_benson_solver.Calculate_solutions(frontier);
    
    std::list<std::pair<std::list<std::string>, Point>> solutions;
    
	std::cout << "Solutions:\n";
    for(auto point : frontier) {
        solutions.push_back(make_pair(std::list<std::string>(), *point));
		std::cout << *point << std::endl;
    }
                            
    add_solutions(solutions.begin(), solutions.end());
    
	std::cout << "ve time: " << dual_benson_solver.vertex_enumeration_time() << "\n";
	std::cout << "lp time: " << solver_time << "\n";
}
    
}

#endif /* defined(__mco__pilp_dual_benson__) */
