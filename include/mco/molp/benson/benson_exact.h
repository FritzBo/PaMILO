#pragma once
/*
 * benson.h
 *
 *  Created on: 20.06.2013
 *      Author: fritz
 */

#ifndef BENSON_H_
#define BENSON_H_

#include <list>

#include <gurobi_c++.h>

#include <mco/molp/basic/molp_model.h>
#include <mco/basic/abstract_solver.h>

namespace mco {
    
class P1;
class P2;
class D2;

template<typename OnlineVertexEnumerator>
class PrimalBensonMolpSolver : public AbstractSolver<std::list<unsigned>> {

public:
	PrimalBensonMolpSolver() = delete;
	PrimalBensonMolpSolver(const mco::MolpModel &model, double epsilon = PrimalBensonMolpSolver::kEpsilon);

	void Init();
	void Solve();

private:

	static constexpr double kEpsilon = 1E-8;

	const double epsilon_;

	Point *p_hat_ = nullptr;

	const mco::MolpModel & model_;
	GRBEnv *grb_env_;
    
    P1 * p1_;
    P2 * p2_;
    D2 * d2_;

	enum State {
		CONSTRUCTED, INITIALIZED, SOLVED
	} state_;

};
    
class SupplementLP {
public:
    virtual GRBModel *model() {
        return grb_model_;
    }
    
protected:
    GRBModel *grb_model_;
    GRBEnv *grb_env_;
    
    const MolpModel &molp_model_;
    
    SupplementLP(GRBEnv *grb_env, const MolpModel &molp_model) : grb_model_(nullptr), grb_env_(grb_env), molp_model_(molp_model) {}
};

class P1 : public SupplementLP {
public:
    P1(GRBEnv *grb_env, const MolpModel &model);
    void set_weights(double * l);
    
private:
    GRBVar * x_;
    GRBConstr * c_;
    double * l_;
    
};

class P2 : public SupplementLP {
public:
    P2(GRBEnv *grb_env, P1 *p1, const MolpModel &model);
    void set_y(Point &y);
    void set_p_hat(Point *p_hat);
    Point * get_x();
    Point * get_u();
    Point * get_l();
    bool is_outsider();
    ~P2();
    
private:
    Point *p_hat_ = nullptr;
    Point * y_;
    GRBVar * x_;
    GRBVar z_;
    GRBConstr * y_constr_;
    GRBConstr * c_;
};

class D2 : public SupplementLP {
public:
    D2(GRBEnv *grb_env, const MolpModel &model);
    void set_y(Point &y);
    Point *get_l();
    Point *get_u();
    
private:
    Point * y_;
    GRBVar *u_;
    GRBVar *l_;
    GRBConstr * x_constr_;
    GRBConstr * z_constr_;
    GRBConstr dual_lock_;
    std::list<GRBConstr> lambda_lock_;
};


} // namespace mco

#endif /* BENSON_H_ */
