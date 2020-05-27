/*
 * ap_benson_dual_solver.cpp
 *
 *  Created on: 14.10.2013
 *      Author: fritz
 */

#include <mco/ap/benson_dual/ap_grb_benson_dual_solver.h>

#include <list>
#include <exception>

using std::list;
using std::exception;

#include <ogdf/basic/Graph.h>

using ogdf::edge;
using ogdf::node;

namespace mco {

template<typename OnlineVertexEnumerator>
APGRBBensonDualSolver<OnlineVertexEnumerator>::APGRBBensonDualSolver(AssignmentInstance & instance,
                                                                     double epsilon)
    :   AbstractAPSolver(instance),
        grb_env_(new GRBEnv()),
        epsilon_(epsilon),
        cycles_(0) {

	grb_env_ = new GRBEnv();
	lp_model_ = new GRBModel(*grb_env_);
	edge e;
	node n;

	variables_ = lp_model_->addVars(instance.graph().numberOfEdges());

	lp_model_->update();

	for(int i = 0; i < instance.graph().numberOfEdges(); ++i) {
		variables_[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
		variables_[i].set(GRB_DoubleAttr_LB, 0);
	}

	forall_nodes(n, instance.graph()) {
		GRBLinExpr c = 0;
		unsigned int i = 0;
		forall_edges(e, instance.graph()) {
			if(e->isIncident(n))
				c += variables_[i];
			++i;
		}
		lp_model_->addConstr(c, GRB_EQUAL, 1, "");
	}

	lp_model_->update();

	dual_benson_solver_ = new DualBensonScalarizer<OnlineVertexEnumerator>(*this, instance.dimension(), epsilon);
}

template<typename OnlineVertexEnumerator>
APGRBBensonDualSolver<OnlineVertexEnumerator>::~APGRBBensonDualSolver() {
	delete lp_model_;
	delete grb_env_;
}

template<typename OnlineVertexEnumerator>
double APGRBBensonDualSolver<OnlineVertexEnumerator>::Solve_scalarization(Point& weights, Point& value) {
	clock_t start = clock();
	edge e;
	GRBLinExpr obj = 0;

	double weighted_value;

	list<GRBConstr> constraints;

	unsigned int i = 0;
	forall_edges(e, *instance().graph()) {
		obj += variables_[i] * (*(*instance().weights())[e] * weights);
		++i;
	}

	lp_model_->setObjective(obj, GRB_MINIMIZE);
	lp_model_->optimize();

	weighted_value = lp_model_->get(GRB_DoubleAttr_ObjVal);
	constraints.push_back(lp_model_->addConstr(obj, GRB_EQUAL, weighted_value, ""));

	for(unsigned int j = 0; j < instance().dimension(); ++j) {
		obj = 0;
		i = 0;
		forall_edges(e, *instance().graph()) {
			obj += variables_[i] * (*(*instance().weights())[e])[j];
			++i;
		}

		lp_model_->setObjective(obj, GRB_MINIMIZE);
		lp_model_->optimize();

		value[j] = lp_model_->get(GRB_DoubleAttr_ObjVal);
		constraints.push_back(lp_model_->addConstr(obj, GRB_EQUAL, value[j], ""));

		lp_model_->update();
	}

	for(auto constr : constraints) {
		try {
			lp_model_->remove(constr);
		} catch(GRBException &e) {
			cout << e.getMessage() << endl;
		}
	}

	cycles_ += clock() - start;

	return weighted_value;
}

} /* namespace mco */
