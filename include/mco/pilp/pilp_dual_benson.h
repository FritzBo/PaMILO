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
#include <typeinfo>
#include <set>
#include <math.h>

#include <mco/pilp/coin.h>

#include <mco/basic/point.h>
#include <mco/basic/abstract_solver.h>
#include <mco/generic/benson_dual/dual_benson_scalarizer.h>

namespace mco {

class ILPSolverAdaptor {
public:
    ILPSolverAdaptor(ILP &ilp, double& time, bool preprocessScale = true) : ilp_(ilp), time_(time) {
		if(!preprocessScale) {
			return;
		}

		int dim = ilp_.dimension;
		
		Point weighting(dim), mini(dim), maxi(dim);
		std::vector<double> spread(dim);
		for(int i = 0; i < dim; i++) {
			weighting[i] = 0;
			mini[i] = std::numeric_limits<double>::max();
			maxi[i] = std::numeric_limits<double>::min();
		}

		for(int i = 0; i < dim; i++) {
			weighting[i] = 1;
			Point value(dim);
			double scalar_value = operator()(weighting, value);
			for(int j = 0; j < dim; j++) {
				mini[j] = std::min(mini[j], value[j]);
				maxi[j] = std::max(maxi[j], value[j]);
			}
			weighting[i] = 0;
		}

		for(int i = 0; i < dim; i++) {
			spread[i] = maxi[i] - mini[i];
		}

		std::nth_element(spread.begin(), spread.begin() + dim/2, spread.end());
		double medianSpread = spread[dim/2];
		for(int i = 0; i < dim; i++) {
			double offset = 0; -maxi[i];
			if(maxi[i] < 0) {
				offset = -maxi[i];
			} else if(mini[i] > 0) {
				offset = -mini[i];
			}
			ilp_.offset[i] = offset;
			ilp_.relScale[i] = exp2(round(log2((maxi[i] - mini[i]) / medianSpread)));
			std::cout << i << ", " << maxi[i] << " " << mini[i] << " spread: " << (maxi[i] - mini[i]) << ", " << log2((maxi[i] - mini[i]) / medianSpread) << std::endl;
		}

		IloNumExprArray objs(ilp_.env);
		IloNumArray weights(ilp_.env);
		IloIntArray prio(ilp_.env);
		IloNumArray absTols(ilp_.env);
		IloNumArray relTols(ilp_.env);

		for(int i = 0; i < ilp_.dimension; i ++) {
			objs.add(ilp_.relScale[i] * (ilp_.obj.getCriterion(i) + ilp_.offset[i]));
			weights.add(1);
			prio.add(0);
			absTols.add(0);
			relTols.add(0);
		}
		auto sense = ilp_.obj.getSense();
		ilp_.model.remove(ilp_.obj);
		ilp_.obj = IloObjective(ilp_.env, IloStaticLex(ilp_.env, objs, weights,
										 prio, absTols, relTols), sense);
		ilp_.model.add(ilp_.obj);
		std::cout << ilp_.obj << std::endl;

		ilp_.cplex.extract(ilp_.model);
	}
    
    inline double operator()(const Point& weighting,
                             Point& value);
    
    inline ~ILPSolverAdaptor();

private:
	ILP &ilp_;

    std::set<Point*, LexPointComparator> known_points_;

	double &time_;
};
    
template<typename OnlineVertexEnumerator>
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
	IloNumExprArray objs(ilp_.env);
	IloNumArray weights(ilp_.env);
	IloIntArray prio(ilp_.env);
	IloNumArray absTols(ilp_.env);
	IloNumArray relTols(ilp_.env);

	for(int i = 0; i < ilp_.dimension; i ++) {
		objs.add(ilp_.obj.getCriterion(i));
		weights.add(weighting[i]);
		prio.add(ilp_.dimension + 1);
		absTols.add(0);
		relTols.add(0);
	}

	for(int i = 0; i < ilp_.dimension; i++) {
		if(weighting[i] < 1E-06) {
			objs.add(ilp_.obj.getCriterion(i));
			weights.add(1);
			prio.add(ilp_.dimension - i);
			absTols.add(0);
			relTols.add(0);
		}
	}

	//std::cout << "objs:\n" << objs << std::endl;
	//std::cout << "weights:\n" << weights << std::endl;
	//std::cout << "prio:\n" << prio << std::endl;
	//std::cout << "absTols:\n" << absTols << std::endl;
	//std::cout << "relTols:\n" << relTols << std::endl;

	auto sense = ilp_.obj.getSense();
	ilp_.model.remove(ilp_.obj);
	ilp_.obj = IloObjective(ilp_.env, IloStaticLex(ilp_.env, objs, weights,
									 prio, absTols, relTols), sense);
	ilp_.model.add(ilp_.obj);
	//std::cout << ilp_.obj << std::endl;

	ilp_.cplex.extract(ilp_.model);
	clock_t start = clock();
	ilp_.cplex.solve();
	time_ +=  (clock() - start) / double(CLOCKS_PER_SEC);

	//TODO: test if bounded

	for(int i = 0; i < ilp_.dimension; i++) {
		value[i] = ilp_.cplex.getValue(objs[i], -1);
	}

    return ilp_.cplex.getMultiObjInfo(IloCplex::MultiObjObjValue, 0);
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
    
	std::cout << "No of solutions:\n";
    for(auto point : frontier) {
        solutions.push_back(make_pair(std::list<std::string>(), *point));
		//std::cout << *point << std::endl;
    }
	std::cout << solutions.size() << std::endl;
                            
    add_solutions(solutions.begin(), solutions.end());
    
	std::cout << "ve time: " << dual_benson_solver.vertex_enumeration_time() << "\n";
	std::cout << "cplex time: " << solver_time << "\n";
}
    
}

#endif /* defined(__mco__pilp_dual_benson__) */
