//
//  pilp_benson_module.cpp
//  pamilo
//
//  Created by Fritz Bökler and Mirko H. Wagner on 28.05.20.
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

        ValueArg<double> epsilon_argument("e", "epsilon", "Epsilon to be used in floating point calculations.", false, 1E-8, "epsilon");

        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");

        cmd.add(epsilon_argument);
        cmd.add(file_name_argument);

        cmd.parse(argc, argv);

        string file_name = file_name_argument.getValue();
        double epsilon = epsilon_argument.getValue();

		ILP ilp;

		LPparser parser;

        parser.getILP(file_name, ilp);

#ifdef USE_CDD
		if(ilp.dimension > 4) {
			PilpDualBensonSolver<OnlineVertexEnumeratorCDD> solver(epsilon);
			solver.Solve(ilp);

			solutions_.insert(solutions_.begin(),
							  solver.solutions().cbegin(),
							  solver.solutions().cend());
		} else
#endif
		{
			PilpDualBensonSolver<GraphlessOVE> solver(epsilon);
			solver.Solve(ilp);

			solutions_.insert(solutions_.begin(),
							  solver.solutions().cbegin(),
							  solver.solutions().cend());
		}


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
