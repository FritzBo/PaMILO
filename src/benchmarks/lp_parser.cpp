//
//  lp_parser.cpp
//  pamilo
//
//  Created by Mirko H. Wagner on 20.06.20
//
//  This file is distributed for academics only
//  under the terms of an MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
//
//

#include <pamilo/benchmarks/lp_parser.h>

#include <string.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

namespace pamilo {
	void LPparser::getILP(string filename, ILP &ilp) {
		try {
			ilp.cplex.importModel(ilp.model, filename.c_str(), ilp.obj, ilp.vars, ilp.cons);
		} catch(IloException &e) {
			cerr << "CPLEX failed to read the file. This is likely, because the input file is corrupted or is not a MOMIP\n";
			exit(-1);
		}
		ilp.cplex.setParam(IloCplex::Param::MultiObjective::Display, 2);
		ilp.cplex.setParam(IloCplex::Param::ParamDisplay, 0);
		ilp.cplex.setOut(ilp.cplexFile);

		ilp.dimension = ilp.obj.getNumCriteria();

		ilp.relScale.resize(ilp.dimension, 1);
		ilp.offset.resize(ilp.dimension, 0);

		ilp.cplex.extract(ilp.model);

		ilp.filename = filename;
	}
}

