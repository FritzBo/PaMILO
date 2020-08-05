//
//  cplex.h
//  pamilo
//
//  Add
//
//  Created by Mirko H. Wagner on 20.06.20
//
//

#pragma once

#include <ilcplex/ilocplex.h>

#include <vector>


class ILP {
public:
	IloEnv env;
	IloModel model;
	IloCplex cplex;
	IloObjective obj;
	IloNumVarArray vars;
	IloRangeArray cons;

	std::vector<double> relScale;
	std::vector<double> offset;

	int dimension;

	ILP() : model(env), cplex(env), vars(env), cons(env), dimension(-1) {}

	~ILP() {
		env.end();
	}
};
