#pragma once
/*
 * ap_benson_dual_solver.h
 *
 *  Created on: 14.10.2013
 *      Author: fritz
 */

#ifndef AP_BENSON_DUAL_SOLVER_H_
#define AP_BENSON_DUAL_SOLVER_H_

#include <functional>

#include <mco/basic/abstract_solver.h>
#include <mco/basic/weight_function_adaptors.h>
#include <mco/ap/basic/abstract_ap_solver.h>
#include <mco/ap/basic/ap_instance.h>
#include <mco/ap/basic/lex_hungarian.h>
#include <mco/generic/benson_dual/dual_benson_scalarizer.h>
#include <mco/generic/benson_dual/ove_fp_v2.h>

namespace mco {

class LexHungarianSolverAdaptor {
public:
    inline LexHungarianSolverAdaptor(AssignmentInstance& ap_instance);
    
    inline double operator()(const Point& weighting, Point& value);
    
private:
    LexHungarianMethod lex_ap_solver_;
    AssignmentInstance& ap_instance_;
};

    
template<typename OnlineVertexEnumerator = GraphlessOVE>
class APBensonDualSolver
: public AbstractSolver<std::list<ogdf::edge>> {
        
public:
    APBensonDualSolver(double epsilon = 1E-8)
    :   epsilon_(epsilon) { }

	void Solve(AssignmentInstance & instance) {
        
		std::list<Point *> frontier;
        
        DualBensonScalarizer<OnlineVertexEnumerator>
        dual_benson_solver_(LexHungarianSolverAdaptor(instance),
                            instance.dimension(),
                            epsilon_);
        
		dual_benson_solver_.Calculate_solutions(frontier);
        
        std::list<std::pair<std::list<ogdf::edge>, Point>> solutions;
        
        for(auto point : frontier) {
            solutions.push_back(make_pair(std::list<ogdf::edge>(), *point));
        }
        
		add_solutions(solutions.begin(), solutions.end());
	}

private:
    double epsilon_;
};
    
        
    
inline LexHungarianSolverAdaptor::
LexHungarianSolverAdaptor(AssignmentInstance& ap_instance)
:   ap_instance_(ap_instance) {
}


inline double LexHungarianSolverAdaptor::
operator()(const Point& weighting, Point& value) {
    
    Point result = lex_ap_solver_.Solve(ap_instance_.graph(),
                                        LexWeightFunctionAdaptor(ap_instance_.graph(),
                                                                 ap_instance_.weights(),
                                                                 weighting),
                                        ap_instance_.dimension() + 1,
                                        ap_instance_.agents());
    
    for(unsigned i = 0; i < ap_instance_.dimension(); ++i) {
        value[i] = result[i + 1];
    }
    
    return result[0];
}




} /* namespace mco */
#endif /* AP_BENSON_DUAL_SOLVER_H_ */
