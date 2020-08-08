//
//  lp_parser.cpp
//  pamilo
//
//  Created by Mirko H. Wagner on 20.06.20
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
		ilp.cplex.importModel(ilp.model, filename.c_str(), ilp.obj, ilp.vars, ilp.cons);
		ilp.cplex.setParam(IloCplex::Param::MultiObjective::Display, 0);
		ilp.cplex.setParam(IloCplex::Param::ParamDisplay, 0);
		ilp.cplex.setOut(ilp.env.getNullStream());

		ilp.dimension = ilp.obj.getNumCriteria();

		ilp.relScale.resize(ilp.dimension, 1);
		ilp.offset.resize(ilp.dimension, 0);

		ilp.cplex.extract(ilp.model);

		ilp.filename = filename;
	}
}
