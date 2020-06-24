//
//  pilp_benson_module.cpp
//  mco
//
//  Created by Fritz Bökler and Mirko H. Wagner on 28.05.20.
//
//

#include "pilp_benson_module.h"

using std::map;
using std::string;
using std::list;
using std::pair;
using std::function;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::ArgException;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;

#include <mco/pilp/pilp_dual_benson.h>
#include <mco/benchmarks/lp_parser.h>
#include <mco/basic/point.h>

using mco::PilpDualBensonSolver;
using mco::Point;
using mco::LPparser;

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
        
        PilpDualBensonSolver<> solver(epsilon);
        
        solver.Solve(ilp);
        
        solutions_.insert(solutions_.begin(),
                          solver.solutions().cbegin(),
                          solver.solutions().cend());
        
    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const list<string>, const Point>>& PilpBensonModule::solutions() {
    return solutions_;
}

string PilpBensonModule::statistics() {
    string stats("");
    return stats;
}