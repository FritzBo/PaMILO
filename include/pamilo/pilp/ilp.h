//
//  ilp.h
//  pamilo
//
//  Created by Mirko H. Wagner on 20.06.20
//
//  This file is distributed under the terms of
//
//  the GNU General Public License v3,
//  a copy of which can be found in the file LICENCE-GPLv3.txt
//
//  OR
//
//  for academics, a MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
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

	std::string filename;
	std::ofstream solFile;
	std::ofstream logFile;

	int dimension;

	ILP() : model(env), cplex(env), vars(env), cons(env), dimension(-1), filename("pilpInstance") {}

	~ILP() {
		env.end();
	}
};
