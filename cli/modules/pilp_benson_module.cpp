//
//  pilp_benson_module.cpp
//  pamilo
//
//  Created by Fritz Bökler and Mirko H. Wagner on 28.05.20.
//
//  This file is distributed for academics only
//  under the terms of an MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
//
//

#include "pilp_benson_module.h"

using std::string;
using std::list;
using std::pair;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::ArgException;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::SwitchArg;

#include <pamilo/pilp/pilp_dual_benson.h>
#include <pamilo/benchmarks/lp_parser.h>
#include <pamilo/basic/point.h>

#include <pamilo/generic/benson_dual/ove_fp_v2.h>

#include <config_autogen.h>
#ifdef USE_CDD
#include <pamilo/generic/benson_dual/ove_cdd.h>
using pamilo::OnlineVertexEnumeratorCDD;
#endif

using pamilo::PilpDualBensonSolver;
using pamilo::GraphlessOVE;
using pamilo::Point;
using pamilo::LPparser;

void PilpBensonModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Dual Benson to find the Pareto-frontier of the parametric integer linear program problem.", ' ', "0.1");

        ValueArg<string> output_name_argument("o", "output", "Basename of the output files. This defaults to <instance>.", false, "", "output");

        ValueArg<double> epsilon_argument("e", "epsilon", "Epsilon to be used in floating point calculations.", false, 1E-7, "epsilon");

        UnlabeledValueArg<string> instance_name_argument("instance", "Name of the instance file.", true, "","instance");

		SwitchArg no_preprocessing_argument("", "no-pre", "Don't run preprocessing. Only use this, if you know all objectives are in rougly the same range and either the lowest or the highest value in each objective is close to 0.", false);

		ValueArg<string> ve_argument("E", "vertex-enumeration", "Which vertex enumeration algorithm is to be used. Options are: cdd, graphless, and auto (default)", false, "auto", "vertex-enumeration");

        cmd.add(output_name_argument);
        cmd.add(epsilon_argument);
        cmd.add(instance_name_argument);
        cmd.add(no_preprocessing_argument);
		cmd.add(ve_argument);

        cmd.parse(argc, argv);

		ILP ilp;

        string instance_name = instance_name_argument.getValue();
        double epsilon = epsilon_argument.getValue();
		string output_name = output_name_argument.getValue();
		bool no_preprocessing = no_preprocessing_argument.getValue();
		string ve = ve_argument.getValue();

		if(output_name == "") {
			output_name == instance_name;
		}
		ilp.solFile.open(output_name + "_sol");
		ilp.logFile.open(output_name + "_log");
		ilp.cplexFile.open(output_name + "_cplex");
		ilp.noPreprocessing = no_preprocessing;

		LPparser parser;

		try {
			parser.getILP(instance_name, ilp);

#ifdef USE_CDD
			if(ve == "cdd" || (ilp.dimension > 4 && ve == "auto")) {
				PilpDualBensonSolver<OnlineVertexEnumeratorCDD> solver(epsilon);
				solver.Solve(ilp);

				solutions_.insert(solutions_.begin(),
								  solver.solutions().cbegin(),
								  solver.solutions().cend());
			} else
#endif
			{
				if(ve == "cdd") {
					std::cerr << "cdd is not activated in cmake!\n";
					exit(0);
				}
				PilpDualBensonSolver<GraphlessOVE> solver(epsilon);
				solver.Solve(ilp);

				solutions_.insert(solutions_.begin(),
								  solver.solutions().cbegin(),
								  solver.solutions().cend());
			}
		} catch (IloException &e) {
			std::cerr << e.getMessage();
		}
		ilp.solFile.close();
		ilp.logFile.close();
		ilp.cplexFile.close();
    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const string, const Point>>& PilpBensonModule::solutions() {
    return solutions_;
}

string PilpBensonModule::statistics() {
    string stats("");
    return stats;
}

