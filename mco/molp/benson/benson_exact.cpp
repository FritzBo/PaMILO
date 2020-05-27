/*
 * benson.cpp
 *
 *  Created on: 20.06.2013
 *      Author: fritz
 */

#include <mco/molp/benson/benson_exact.h>

#include <iostream>
#include <exception>
#include <set>
#include <list>
#include <cmath>
#include <cassert>

using std::invalid_argument;
using std::exception;
using std::set;
using std::list;
using std::abs;

#include <gurobi_c++.h>
#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::randomDoubleNormal;
using ogdf::randomDouble;

#include <mco/molp/basic/molp_model.h>

namespace mco {

template<typename OnlineVertexEnumerator>
PrimalBensonMolpSolver<OnlineVertexEnumerator>::PrimalBensonMolpSolver(const MolpModel &model, double epsilon) :
		epsilon_(epsilon), p_hat_(nullptr), model_(model), grb_env_(nullptr), p1_(nullptr), p2_(nullptr), d2_(nullptr) {

	state_ = CONSTRUCTED;
}

template<typename OnlineVertexEnumerator>
void PrimalBensonMolpSolver<OnlineVertexEnumerator>::Init() {
	if(state_ != CONSTRUCTED)
		return;

	grb_env_ = new GRBEnv();

	p1_ = new P1(grb_env_, model_);
	p2_ = new P2(grb_env_, p1_, model_);
	d2_ = new D2(grb_env_, model_);

	state_ = INITIALIZED;
}

template<typename OnlineVertexEnumerator>
void PrimalBensonMolpSolver<OnlineVertexEnumerator>::Solve() {
	if(state_ == CONSTRUCTED)
		return;

	unsigned int dimension = model_.objectives();
	unsigned int variables = model_.variables();
	Point b = Point(model_.b(), model_.constraints());

	/*
	 * Finding ideal point
	 */
	double *s0 = new double[dimension];
	double *l = new double[dimension];
	for(unsigned int i = 0; i < dimension; ++i)
		l[i] = i == 0 ? 1 : 0;

	for(unsigned int i = 0; i < dimension; ++i) {
		p1_->set_weights(l);
		p1_->model()->optimize();
		if(p1_->model()->get(GRB_IntAttr_Status) != GRB_OPTIMAL)
			throw new invalid_argument("Unbounded or infeasible model.");

		s0[i] = p1_->model()->get(GRB_DoubleAttr_ObjVal);

		if(i < model_.objectives() - 1) {
			l[i] = 0;
			l[i + 1] = 1;
		}
	}

	Point ideal_point(s0, dimension);

	std::cout << "Ideal point is: " << ideal_point << std::endl;

	delete[] s0;
	delete[] l;

	// Done finding ideal point

	OnlineVertexEnumerator vertex_container(ideal_point, dimension, epsilon_);

	unsigned int iterations = 0;

	while(vertex_container.has_next()) {
		iterations += 1;

		Point *s = vertex_container.next_vertex();

		std::cout << "******************************************************************************************" << std::endl;
		std::cout << "******************************************************************************************" << std::endl;

		cout << "Iteration: " << iterations << endl;

		std::cout << "s is: " << *s << std::endl;

		/*
		 * Check if s in Y
		 */
		Point *x, *y, *y_prime;

		if(p_hat_ == nullptr) {

			p_hat_ = new Point(dimension);
			for(unsigned int i = 0; i < dimension; ++i)
				(*p_hat_)[i] = (*s)[i] + randomDouble(1, 5);

			p2_->set_p_hat(p_hat_);

			p2_->set_y(*s);
			p2_->model()->write("test.lp");
			p2_->model()->optimize();

			x = p2_->get_x();

			y = new Point(dimension);
			for(unsigned int i = 0; i < dimension; ++i) {
				for(unsigned int j = 0; j < variables; ++j)
					(*y)[i] += (*x)[j] * model_.C()[model_.variables() * i + j];
			}

			y_prime = y;

			Point perturbation(dimension);
			for(unsigned int i = 0; i < dimension; ++i)
				perturbation[i] = randomDouble(0.8, 3);

			p_hat_ = new Point((*y) + perturbation * p2_->model()->get(GRB_DoubleAttr_ObjVal) * 0.5);
			cout << "p_hat: " << *p_hat_ << endl;
			p2_->set_p_hat(p_hat_);

		} else {

			p2_->set_y(*s);
			p2_->model()->write("test.lp");
			p2_->model()->optimize();

			x = p2_->get_x();

			y_prime = new Point(dimension);
			for(unsigned int i = 0; i < dimension; ++i) {
				for(unsigned int j = 0; j < variables; ++j)
					(*y_prime)[i] += (*x)[j] * model_.C()[model_.variables() * i + j];
			}

			double * y_values = new double[dimension];
			for(unsigned int i = 0; i < dimension; ++i)
				y_values[i] = (*s)[i] + p2_->model()->get(GRB_DoubleAttr_ObjVal) * ((*p_hat_)[i] - (*s)[i]);

			y = new Point(y_values, dimension);
			delete[] y_values;
		}

		// We only deal with points OUTSIDE the outcome polyhedron
		assert(p2_->model()->get(GRB_DoubleAttr_ObjVal) > -epsilon_);

		// Vertex is inside S and thus an efficient vertex, proceed
		if(p2_->model()->get(GRB_DoubleAttr_ObjVal) < epsilon_) {
			add_solution(y);
			continue;
		}

		// Finished checking

		std::cout << "x is: " << (*x) << std::endl;
		std::cout << "y is: " << (*y) << std::endl;
		std::cout << "y' is: " << (*y_prime) << std::endl;

		if(p2_->is_outsider())
			std::cout << "s is an outsider." << std::endl;

		/*
		 * Find lexicographic minimal lambda
		 */
//		d2_->set_y(*y);
//
//		for(unsigned int i = 0; i < dimension - 1; ++i) {
//		}

		// Finished finding lexicographic minimal lambda

		Point *l = p2_->get_l();
		Point *u = p2_->get_u();

		cout << "l is: " << *l << endl;
		cout << "u is: " << *u << endl;

		vertex_container.add_hyperplane(*s, *l, b * (*u));

		delete s;
		delete x;
		delete y;
		delete l;
		delete u;
	}

	cout << "Iteration count: " << iterations << endl;
	cout << "Number of inequalities: " << vertex_container.number_of_hyperplanes() << endl;
}

P1::P1(GRBEnv *grb_env, const MolpModel &model) :
				SupplementLP(grb_env, model), l_(nullptr) {

	grb_model_ = new GRBModel(*grb_env);
	x_ = grb_model_->addVars(molp_model_.variables());
	grb_model_->update();
	for(unsigned int i = 0; i < model.variables(); ++i) {
		x_[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
		x_[i].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
	}

	c_ = new GRBConstr[molp_model_.constraints()];
	for(unsigned int i = 0; i < model.constraints(); ++i) {
		GRBLinExpr lhs = 0;
		for(unsigned int j = 0; j < model.variables(); ++j) {
			lhs += x_[j] * molp_model_.A()[molp_model_.variables() * i + j];
		}
		c_[i] = grb_model_->addConstr(lhs, GRB_GREATER_EQUAL, molp_model_.b()[i]);
	}

	grb_model_->update();
	grb_model_->write("test.lp");
}

void P1::set_weights(double *l) {
	l_ = l;

	GRBLinExpr obj = 0;
	for(unsigned int i = 0; i < molp_model_.objectives(); ++i) {
		for(unsigned int j = 0; j < molp_model_.variables(); ++j) {
			obj += l_[i] * molp_model_.C()[molp_model_.variables() * i + j] * x_[j];
		}
	}
	grb_model_->setObjective(obj, GRB_MINIMIZE);
}

P2::P2(GRBEnv *grb_env, P1 *p1, const MolpModel &molp_model) :
		SupplementLP(grb_env, molp_model), p_hat_(nullptr) {

	y_ = nullptr;
	y_constr_ = nullptr;
	grb_model_ = new GRBModel(*p1->model());
	x_ = grb_model_->getVars();
	c_ = grb_model_->getConstrs();
	z_ = grb_model_->addVar(-GRB_INFINITY, GRB_INFINITY, 1, GRB_CONTINUOUS);
	grb_model_->update();

	grb_model_->setObjective(z_ + 0, GRB_MINIMIZE);
}

void P2::set_p_hat(Point *p_hat) {
	p_hat_ = p_hat;
}

void P2::set_y(Point &y) {
	if(p_hat_ == nullptr)
		return;

	y_ = &y;

	if(y_constr_ != nullptr)
		for(unsigned int i = 0; i < molp_model_.objectives(); ++i)
			grb_model_->remove(y_constr_[i]);

	y_constr_ = new GRBConstr[molp_model_.objectives()];

	double direction_norm = 0;
	for(unsigned int i = 0; i < molp_model_.objectives(); ++i)
		direction_norm += abs((*p_hat_)[i] - (*y_)[i]);

	for(unsigned int i = 0; i < molp_model_.objectives(); ++i) {
		GRBLinExpr lhs = 0;
		for(unsigned int j = 0; j < molp_model_.variables(); ++j)
			lhs += molp_model_.C()[molp_model_.variables() * i + j] * x_[j];

		double normed_direction = 1 * ((*p_hat_)[i] - (*y_)[i]) / direction_norm;
		y_constr_[i] = grb_model_->addConstr(lhs, GRB_LESS_EQUAL, (*y_)[i] + z_ * normed_direction);
//		y_constr_[i] = grb_model_->addConstr(lhs, GRB_LESS_EQUAL, (*y_)[i] + z_ * ((*p_hat_)[i] - (*y_)[i]));
	}
	grb_model_->update();
}

Point * P2::get_x() {
	double *x = new double[molp_model_.variables()];

	for(unsigned int i = 0; i < molp_model_.variables(); ++i)
		x[i] = x_[i].get(GRB_DoubleAttr_X);

	Point * point = new Point(x, molp_model_.variables());
	delete[] x;

	return point;
}

Point * P2::get_l() {
	double *l = new double[molp_model_.objectives()];

	for(unsigned int i = 0; i < molp_model_.objectives(); ++i)
		l[i] = -y_constr_[i].get(GRB_DoubleAttr_Pi);

	Point * point = new Point(l, molp_model_.objectives());
	delete[] l;

	return point;
}

Point * P2::get_u() {
	double *u = new double[molp_model_.constraints()];

	for(unsigned int i = 0; i < molp_model_.constraints(); ++i)
		u[i] = c_[i].get(GRB_DoubleAttr_Pi);

	Point * point = new Point(u, molp_model_.constraints());
	delete[] u;

	return point;
}

bool P2::is_outsider() {
	for(unsigned int i = 0; i < molp_model_.objectives(); ++i)
		if(y_constr_[i].get(GRB_DoubleAttr_Slack) > 1E-6)
			return true;
	return false;
}

P2::~P2() {
	cout << "p2 destruct" << endl;
	delete p_hat_;
}

D2::D2(GRBEnv *grb_env, const MolpModel &model)
:   SupplementLP(grb_env, model),
    y_(nullptr),
    z_constr_(nullptr),
    lambda_lock_() {

	grb_model_ = new GRBModel(*grb_env_);

	u_ = grb_model_->addVars(molp_model_.constraints());
	grb_model_->update();

	for(unsigned int i = 0; i < molp_model_.constraints(); ++i) {
		u_[i].set(GRB_DoubleAttr_LB, 0);
		u_[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
	}

	l_ = grb_model_->addVars(molp_model_.objectives());
	grb_model_->update();
	for(unsigned int i = 0; i < molp_model_.objectives(); ++i) {
		l_[i].set(GRB_DoubleAttr_LB, 0);
		l_[i].set(GRB_DoubleAttr_UB, GRB_INFINITY);
	}

	x_constr_ = new GRBConstr[molp_model_.variables()];
	for(unsigned int i = 0; i < molp_model_.variables(); ++i) {
		GRBLinExpr lhs = 0;
		for(unsigned int j = 0; j < molp_model_.constraints(); ++j)
			lhs += molp_model_.A()[molp_model_.variables() * j + i] * u_[j];
		for(unsigned int j = 0; j < molp_model_.objectives(); j++)
			lhs -= molp_model_.C()[molp_model_.variables() * j + i] * l_[j];
		x_constr_[i] = grb_model_->addConstr(lhs, GRB_EQUAL, 0);
	}

	GRBLinExpr l_sum = 0;
	for(unsigned int i = 0; i < molp_model_.objectives(); ++i)
		l_sum += l_[i];
	grb_model_->addConstr(l_sum, GRB_EQUAL, 1);
}

void D2::set_y(Point &y) {
	y_ = &y;

	GRBLinExpr obj = 0;
	for(unsigned int i = 0; i < molp_model_.constraints(); ++i)
		obj += molp_model_.b()[i] * u_[i];
	for(unsigned int i = 0; i < molp_model_.objectives(); ++i)
		obj -= (*y_)[i] * l_[i];
	dual_lock_ = grb_model_->addConstr(obj, GRB_GREATER_EQUAL, 0);
	grb_model_->update();
}

Point * D2::get_l() {
	double *l = new double[molp_model_.objectives()];

	for(unsigned int i = 0; i < molp_model_.objectives(); ++i)
		l[i] = l_[i].get(GRB_DoubleAttr_X);

	Point * point = new Point(l, molp_model_.objectives());
	delete[] l;

	return point;

}

Point * D2::get_u() {
	double *u = new double[molp_model_.constraints()];

	for(unsigned int i = 0; i < molp_model_.constraints(); ++i)
		u[i] = u_[i].get(GRB_DoubleAttr_X);

	Point * point = new Point(u, molp_model_.constraints());
	delete[] u;

	return point;
}

} // namespace mco
