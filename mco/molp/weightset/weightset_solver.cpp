/*
 * weightset_solver.cpp
 *
 *  Created on: 19.07.2013
 *      Author: fritz
 */

#include <list>
#include <queue>
#include <set>
#include <cmath>
#include <utility>

using std::list;
using std::deque;
using std::set;
using std::abs;
using std::pair;
using std::make_pair;

#include <gurobi_c++.h>
#include <setoper.h>
#include <cdd.h>

#include <mco/molp/weightset/weightset_solver.h>
#include <mco/molp/basic/molp_model.h>

namespace mco {

WeightSetMolpSolver::WeightSetMolpSolver(const MolpModel& model) : molp_model_(model), grb_env_(new GRBEnv()), grb_model_(nullptr), x_(nullptr), c_(nullptr) {
}

WeightSetMolpSolver::WeightSetMolpSolver(GRBEnv* grb_env, const MolpModel& model) : molp_model_(model), grb_env_(grb_env), grb_model_(nullptr), x_(nullptr), c_(nullptr) {
}

void WeightSetMolpSolver::Init() {
	grb_model_ = new GRBModel(*grb_env_);

	x_ = grb_model_->addVars(molp_model_.variables());
	grb_model_->update();
	for(unsigned int i = 0; i < molp_model_.variables(); ++i) {
		x_[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
		x_[i].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
	}

	c_ = new GRBConstr[molp_model_.constraints()];
	for(unsigned int i = 0; i < molp_model_.constraints(); ++i) {
		GRBLinExpr lhs = 0;
		for(unsigned int j = 0; j < molp_model_.variables(); ++j) {
			lhs += x_[j] * molp_model_.A()[molp_model_.variables() * i + j];
		}
		c_[i] = grb_model_->addConstr(lhs, GRB_GREATER_EQUAL, molp_model_.b()[i]);
	}

	grb_model_->update();
}

class LexComparator {
	unsigned int dimension_;
public:
	LexComparator(unsigned int dimension) : dimension_(dimension) {}

	bool operator()(double * label1, double * label2) {
		for(unsigned int i = 0; i < dimension_; ++i) {
			if(label1[i] - label2[i] < -1E-6) {
				return true;
			} else if(label1[i] - label2[i] > 1E-6) {
				return false;
			}
		}
		return false;
	}
};

void WeightSetMolpSolver::Solve() {
	unsigned int objectives = molp_model_.objectives();
	unsigned int variables = molp_model_.variables();
	unsigned int constraints = molp_model_.constraints();

	double *l, *x, *y, *constr, *obj_constr;

	Init();

	dd_set_global_constants();

	/* ************************************************
	 * Finding an initial solution/extreme point pair *
	 *************************************************/
	l = new double[objectives];
	for(unsigned int i = 0; i < objectives; ++i)
		l[i] = 1;

	x = new double[variables];
	y = new double[objectives];
	lexopt_on_face(l, GRB_MINIMIZE, variables, objectives, x, y);

	delete[] l;

	std::cout << "x is: (" << x[0];
	for(unsigned int i = 1; i < variables; ++i)
		std::cout << ", " << x[i];
	std::cout << ")" << std::endl;

	std::cout << "y is: (" << y[0];
	for(unsigned int i = 1; i < objectives; ++i)
		std::cout << ", " << y[i];
	std::cout << ")" << std::endl;

	// Done finding initial solution/extreme point pair

	deque<pair<double *, double *>> pairs;

	LexComparator comp(objectives);
	set<double *, std::function<bool (double *, double *)>> permanent_points(comp);

	auto initial_pair = make_pair(x,y);
	pairs.emplace_back(initial_pair);

	permanent_points.insert(y);

	while(!pairs.empty()) {

		x = pairs.front().first;
		y = pairs.front().second;

		pairs.pop_front();

		/* ***********************************************************
		 * Finding all equality constraints and applying cost matrix *
		 ************************************************************/
		list<double *> tangent_cone_constraints;

		for(unsigned int i = 0; i < constraints; ++i) {

			constr = new double[variables];

			double lh_value = 0;
			for(unsigned int j = 0; j < variables; ++j) {
				constr[j] = molp_model_.A()[variables * i + j];
				lh_value += x[j] * constr[j];
			}

			if(abs(lh_value - molp_model_.b()[i]) < 1E-6) {
				obj_constr = new double[objectives];

				for(unsigned int j = 0; j < objectives; ++j) {
					obj_constr[j] = 0;

					for(unsigned int k = 0; k < variables; ++k)
						obj_constr[j] += constr[k] * molp_model_.C()[variables * j + k];
				}

				tangent_cone_constraints.push_back(obj_constr);

				delete[] constr;
			}
		}

		assert(tangent_cone_constraints.size() >= variables);

		for(auto constr: tangent_cone_constraints) {
			std::cout << "(" << constr[0];
			for(unsigned int i = 1; i < objectives; ++i)
				std::cout << ", " << constr[i];
			std::cout << ")" << std::endl;
		}

		// Done finding intial inequalities in objective space

		/* *******************************************
		 * Finding interior point of the normal cone *
		 *********************************************/

		double *normal_vector = new double[objectives];

		cone_combination(tangent_cone_constraints, normal_vector, objectives);

		std::cout << "normal vector: (" << normal_vector[0];
		for(unsigned int i = 1; i < objectives; ++i)
			std::cout << ", " << normal_vector[i];
		std::cout << ")" << std::endl;

		// Done finding normal vector

		/* ********************************************
		 * Constructing points on edges adjacent to y *
		 *********************************************/

		dd_ErrorType err;
		dd_MatrixPtr S = dd_CreateMatrix(tangent_cone_constraints.size() + 1, objectives + 1);

		unsigned int i = 0;
		for(auto constr: tangent_cone_constraints) {
			double value = 0;
			for(unsigned int j = 0; j < objectives; ++j) {
				value += y[j] * constr[j];
				dd_set_d(S->matrix[i][j + 1], constr[j]);
			}
			dd_set_d(S->matrix[i][0], -value);

			i++;
		}

		double normal_value = 0;
		for(unsigned int i = 0; i < objectives; ++i) {
			normal_value += normal_vector[i] * y[i];
			dd_set_d(S->matrix[tangent_cone_constraints.size()][i + 1], normal_vector[i]);
		}
		dd_set_d(S->matrix[tangent_cone_constraints.size()][0], -(normal_value + 1));

		S->representation = dd_Inequality;

		dd_PolyhedraPtr P = dd_DDMatrix2Poly(S, &err);
		dd_MatrixPtr G = dd_CopyGenerators(P);

		dd_WriteMatrix(stdout, G);
		dd_WriteMatrix(stdout, S);

		list<double *> edge_points;

		for(int i = 0; i < G->rowsize; ++i) {
			if(dd_get_d(G->matrix[i][0]) != 1)
				continue;

			obj_constr = new double[objectives];
			for(unsigned int j = 0; j < objectives; ++j)
				obj_constr[j] = dd_get_d(G->matrix[i][j + 1]);

			if(is_dominated(obj_constr, y, objectives, 1E-6)) {
				delete[] obj_constr;
				continue;
			}

			edge_points.push_back(obj_constr);
		}

		for(auto point: edge_points) {
			std::cout << "(" << point[0];
			for(unsigned int i = 1; i < objectives; ++i)
				std::cout << ", " << point[i];
			std::cout << ")" << std::endl;
		}

		// Done finding edge points

		/* ****************************
		 * Finding new extreme points *
		 *****************************/

		double *point;
		while(!edge_points.empty()) {
			point = edge_points.front();
			edge_points.pop_front();

			list<double *> equality_constr;

			for(double *tangent_constr: tangent_cone_constraints)
				if(abs(inner_product(tangent_constr, point, objectives) - inner_product(tangent_constr, y, objectives)) < 1E-6)
					equality_constr.push_back(tangent_constr);

			cone_combination(equality_constr, normal_vector, objectives);

			std::cout << "normal: (" << normal_vector[0];
			for(unsigned int i = 1; i < objectives; ++i)
				std::cout << ", " << normal_vector[i];
			std::cout << ")" << std::endl;

			double *new_x = new double[variables];
			double *new_y = new double[objectives];

			lexopt_on_face(normal_vector, GRB_MAXIMIZE, variables, objectives, new_x, new_y);

			if(is_dominated(new_y, y, objectives, 1E-6)) {
				std::cout << "is dominated" << std::endl;
				delete[] new_x;
				delete[] new_y;
				continue;
			}

			if(!epsilon_equal(new_y, y, objectives, 1E-6)) {
				std::cout << "Adding new point: (" << new_y[0];
				for(unsigned int i = 1; i < objectives; ++i)
					std::cout << ", " << new_y[i];
				std::cout << ")" << std::endl;

				if(permanent_points.count(new_y) == 0) {
					pairs.emplace_back(make_pair(new_x, new_y));
					permanent_points.insert(new_y);
				} else {
					delete[] new_x;
					delete[] new_y;
				}
			} else {
				lexopt_on_face(normal_vector, GRB_MINIMIZE, variables, objectives, new_x, new_y);

				if(is_dominated(new_y, y, objectives, 1E-6))
					continue;

				if(permanent_points.count(new_y) == 0) {
					pairs.emplace_back(make_pair(new_x, new_y));
					permanent_points.insert(new_y);
				} else {
					delete[] new_x;
					delete[] new_y;
				}
			}

		}

	}

	dd_free_global_constants();

	std::cout << "Finished. Extreme points are:" << std::endl;
	for(auto point: permanent_points) {
		std::cout << "(" << point[0];
		for(unsigned int i = 1; i < objectives; ++i)
			std::cout << ", " << point[i];
		std::cout << ")" << std::endl;

	}
}

double WeightSetMolpSolver::inner_product(double* vector1, double* vector2,
		unsigned int dimension) {
	double result = 0;

	for(unsigned int i = 0; i < dimension; ++i)
		result += vector1[i] * vector2[i];

	return result;
}

void WeightSetMolpSolver::cone_combination(list<double *> vertices, double *combination, unsigned int dimension) {
	for(auto vertex: vertices) {
		for(unsigned int j = 0; j < dimension; ++j)
			if(vertex == *vertices.begin())
				//combination[j] = 1 / (double) vertices.size() * vertex[j];
				combination[j] = vertex[j];
			else
				//combination[j] += 1 / (double) vertices.size() * vertex[j];
				combination[j] += vertex[j];
	}
}

void WeightSetMolpSolver::set_weights(double *l) {
	GRBLinExpr obj = 0;
	for(unsigned int i = 0; i < molp_model_.objectives(); ++i) {
		for(unsigned int j = 0; j < molp_model_.variables(); ++j) {
			obj += l[i] * molp_model_.C()[molp_model_.variables() * i + j] * x_[j];
		}
	}
	grb_model_->setObjective(obj, GRB_MINIMIZE);
}

void WeightSetMolpSolver::lock_weight(double* l, double value) {
	GRBLinExpr lhs = 0;
	for(unsigned int i = 0; i < molp_model_.objectives(); ++i) {
		for(unsigned int j = 0; j < molp_model_.variables(); ++j) {
			lhs += l[i] * molp_model_.C()[molp_model_.variables() * i + j] * x_[j];
		}
	}

	weight_lock_ = grb_model_->addConstr(lhs, GRB_EQUAL, value);
}

void WeightSetMolpSolver::lock_objective_value(unsigned int objective, double value) {
	GRBLinExpr lhs = 0;

	for(unsigned int i = 0; i < molp_model_.variables(); ++i)
		lhs += x_[i] * molp_model_.C()[molp_model_.variables() * objective + i];

	lock_constraints_.emplace_back(grb_model_->addConstr(lhs, GRB_EQUAL, value));
}

void WeightSetMolpSolver::optimize_objective(unsigned int objective, char sense) {
	GRBLinExpr obj = 0;
	for(unsigned int i = 0; i < molp_model_.variables(); ++i)
		obj += x_[i] * molp_model_.C()[molp_model_.variables() * objective + i];
	grb_model_->setObjective(obj, sense);
}

void WeightSetMolpSolver::unlock() {
	grb_model_->remove(weight_lock_);
	for(auto constr: lock_constraints_)
		grb_model_->remove(constr);

	lock_constraints_.clear();
}

// TODO: Aufteilen: lexopt und lock on face
void WeightSetMolpSolver::lexopt_on_face(double *l, char sense, unsigned int variables, unsigned int objectives, double *x, double *y) {
	set_weights(l);
	grb_model_->optimize();

	assert(grb_model_->get(GRB_IntAttr_Status) == GRB_OPTIMAL);

	lock_weight(l, grb_model_->get(GRB_DoubleAttr_ObjVal));

	optimize_objective(0, sense);
	grb_model_->optimize();

	assert(grb_model_->get(GRB_IntAttr_Status) == GRB_OPTIMAL);

	y[0] = grb_model_->get(GRB_DoubleAttr_ObjVal);


	for(unsigned int i = 1; i < objectives; ++i) {
		lock_objective_value(i - 1, y[i - 1]);

		optimize_objective(i, sense);
		grb_model_->optimize();

		assert(grb_model_->get(GRB_IntAttr_Status) == GRB_OPTIMAL);

		y[i] = grb_model_->get(GRB_DoubleAttr_ObjVal);
	}

	for(unsigned int i = 0; i < variables; ++i) {
		x[i] = x_[i].get(GRB_DoubleAttr_X);
	}

	unlock();
}

bool WeightSetMolpSolver::epsilon_equal(double* point1, double* point2,
		unsigned int dimension, double epsilon = 1E-6) {

	for(unsigned int i = 0; i < dimension; ++i)
		if(abs(point1[i] - point2[i]) > epsilon)
			return false;

	return true;
}

bool WeightSetMolpSolver::is_dominated(double* point1, double* point2,
		unsigned int dimension, double epsilon = 1E-6) {

	for(unsigned int i = 0; i < dimension; ++i)
		if(point1[i] - point2[i] < epsilon)
			return false;

	return true;
}

}
