//
//  pilp_dual_benson.h
//  pamilo
//
//  Created by Fritz Bökler and Mirko H. Wagner on 28.05.20.
//
//  This file is distributed for academics only
//  under the terms of an MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
//
//

#pragma once

#include <functional>
#include <typeinfo>
#include <set>
#include <math.h>
#include <algorithm>
#include <limits>

#include <pamilo/pilp/ilp.h>

#include <pamilo/basic/point.h>
#include <pamilo/basic/abstract_solver.h>
#include <pamilo/generic/benson_dual/dual_benson_scalarizer.h>

namespace pamilo {

struct Log {
	double solver_time = 0;
	int nSolves = 0;
};

class ILPSolverPrinter {
public:
	ILPSolverPrinter(ILP &ilp) : ilp_(ilp) {}
    inline void operator()(std::pair<std::string, Point*>);

private:
	ILP &ilp_;
};

void ILPSolverPrinter::operator()(std::pair<std::string, Point*> sol) {
	Point &point = *(sol.second);
	Point pointScaled(ilp_.dimension);
	for(int i = 0; i < ilp_.dimension; i++) {
		pointScaled[i] = (point[i] / ilp_.relScale[i]) - ilp_.offset[i];
	}
	ilp_.solFile << "[ " << pointScaled << " ]" << sol.first;
#ifdef TEST
	std::cout << "1 ";
#endif
	std::cout << pointScaled << std::endl;
}

class ILPSolverAdaptor {
public:
    ILPSolverAdaptor(ILP &ilp, Log& log)
		: ilp_(ilp), log_(log)
	{
		if(ilp.noPreprocessing) {
			ilp.logFile << "preprocessing time: 0\n";
			return;
		}

		clock_t start = clock();
		int dim = ilp_.dimension;

		Point weighting(dim), mini(dim), maxi(dim);
		std::vector<double> spread(dim);
		for(int i = 0; i < dim; i++) {
			weighting[i] = 0;
			mini[i] = std::numeric_limits<double>::max();
			maxi[i] = std::numeric_limits<double>::lowest();
		}

		for(int i = 0; i < dim; i++) {
			weighting[i] = 1;
			Point value(dim);
			std::string sol;
			double scalar_value = operator()(weighting, value, sol);
			for(int j = 0; j < dim; j++) {
				mini[j] = std::min(mini[j], value[j]);
				maxi[j] = std::max(maxi[j], value[j]);
			}
			weighting[i] = 0;
		}

		int zeroCnt = 0;
		for(int i = 0; i < dim; i++) {
			if(maxi[i] == std::numeric_limits<double>::max()
					|| mini[i] == std::numeric_limits<double>::lowest())
			{
				std::cout << "not bounded in objective " << i << std::endl;
				exit(0);
			}
			spread[i] = maxi[i] - mini[i];
			zeroCnt += (spread[i] == 0 ? 1 : 0);
		}

		auto sense = ilp_.obj.getSense();

		std::nth_element(spread.begin(), spread.begin() + (dim+zeroCnt)/2, spread.end());
		double medianSpread = spread[(dim+zeroCnt)/2];
		for(int i = 0; i < dim; i++) {
			double offset = (maxi[i] < 0 ? -maxi[i] : -mini[i]);
			ilp_.offset[i] = offset;
			if(maxi[i] - mini[i] > 0 || medianSpread == 0) {
				ilp_.relScale[i] = exp2(round(log2(medianSpread / (maxi[i] - mini[i]))));
			}
			ilp_.relScale[i] *= sense;
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
			relTols.add(2); }
		ilp_.model.remove(ilp_.obj);
		ilp_.obj = IloObjective(ilp_.env, IloStaticLex(ilp_.env, objs, weights,
										 prio, absTols, relTols), IloObjective::Minimize);
		ilp_.model.add(ilp_.obj);

		ilp_.cplex.extract(ilp_.model);
		ilp.logFile << "preprocessing time: " << (clock() - start) / (double) CLOCKS_PER_SEC << "\n";
	}

    inline double operator()(const Point& weighting,
                             Point& value,
							 std::string &sol);

    inline ~ILPSolverAdaptor();

private:
	ILP &ilp_;

    std::set<Point*, LexPointComparator> known_points_;

	Log &log_;
};

template<typename OnlineVertexEnumerator>
class PilpDualBensonSolver : public AbstractSolver<std::string> {
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
           Point& value,
		   std::string &sol)
{
	ilp_.logFile << "weighting " << weighting << std::endl;
	IloNumExprArray objs(ilp_.env);
	IloNumArray weights(ilp_.env);
	IloIntArray prio(ilp_.env);
	IloNumArray absTols(ilp_.env);
	IloNumArray relTols(ilp_.env);

	auto sense = ilp_.obj.getSense();

	for(int i = 0; i < ilp_.dimension; i ++) {
		objs.add(ilp_.obj.getCriterion(i) * sense);
		weights.add(weighting[i]);
		prio.add(2);
		absTols.add(0);
		relTols.add(0);
	}

	for(int i = 0; i < ilp_.dimension; i++) {
		if(weighting[i] < 1E-06) {
			objs.add(ilp_.obj.getCriterion(i) * sense);
			weights.add(1);
			prio.add(1);
			absTols.add(0);
			relTols.add(0);
		}
	}

	ilp_.model.remove(ilp_.obj);
	ilp_.obj = IloObjective(ilp_.env, IloStaticLex(ilp_.env, objs, weights,
									 prio, absTols, relTols), IloObjective::Minimize);
	ilp_.model.add(ilp_.obj);

	ilp_.cplex.extract(ilp_.model);
	bool doRun = true;
	clock_t start = clock();

	ilp_.cplex.solve();
	log_.solver_time +=  (clock() - start) / double(CLOCKS_PER_SEC);
	log_.nSolves++;

	//TODO: test if bounded (is this still current?)

	for(int i = 0; i < ilp_.dimension; i++) {
		value[i] = ilp_.cplex.getValue(objs[i], -1);
	}

	for(int i = 0; i < ilp_.vars.getSize(); i++) {
		auto var = ilp_.vars[i];
		sol += " ";
		sol += var.getName();
		sol += "=";
		if(var.getType() != 2) {
			sol += std::to_string(int(ilp_.cplex.getValue(var)));
		} else {
			sol += std::to_string(ilp_.cplex.getValue(var));
		}
	}
	sol += "\n";

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
	clock_t start = clock();
	ilp.startTime = start;

    std::list<std::pair<std::string, Point *>> frontier;
	Log log;

    DualBensonScalarizer<OnlineVertexEnumerator, std::string>
    dual_benson_solver(ILPSolverAdaptor(ilp, log),
                       ILPSolverPrinter(ilp),
                       ilp.dimension,
                       epsilon_);

    dual_benson_solver.Calculate_solutions(frontier);

//#define TEST
#ifdef TEST
	for(int i = 0; i < ilp.dimension; i++) {
		std::cout << "0";
		for(int j = 0; j < ilp.dimension; j++) {
			if(i != j) {
				std::cout << " 0";
			} else {
				std::cout << " 1";
			}
		}
		std::cout << std::endl;
	}
#endif

    for(auto sol : frontier) {
		Point &point = *(sol.second);
		Point pointScaled(ilp.dimension);
		for(int i = 0; i < ilp.dimension; i++) {
			pointScaled[i] = (point[i] / ilp.relScale[i]) - ilp.offset[i];
		}
		delete sol.second;
		add_solution(sol.first, pointScaled);
    }

	ilp.logFile << "time: " << (clock() - start) / (double) CLOCKS_PER_SEC << "\n";
	ilp.logFile << "vertex enumeration time: " << dual_benson_solver.vertex_enumeration_time() << "\n";
	ilp.logFile << "cplex time: " << log.solver_time << "\n";
	ilp.logFile << "cplex solves: " << log.nSolves << "\n";
	ilp.logFile << "time per solve: " << log.solver_time/log.nSolves << "\n";
}
}

