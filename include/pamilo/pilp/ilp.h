//
//  ilp.h
//  pamilo
//
//  Created by Mirko H. Wagner on 20.06.20
//
//  This file is distributed for academics only
//  under the terms of an MIT license based license,
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
	std::ofstream cplexFile;
	bool noPreprocessing;

	int startTime;

	int dimension;

	ILP() :
		model(env),
		cplex(env),
		vars(env),
		cons(env),
		dimension(-1),
		filename("pilpInstance"),
		solFile(""),
		logFile(""),
		cplexFile(""),
		noPreprocessing(false),
		startTime(-1)
	{}

	~ILP() {
		env.end();
	}
};

