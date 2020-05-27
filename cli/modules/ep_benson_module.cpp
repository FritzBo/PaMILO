//
//  benson_module.cpp
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#include "ep_benson_module.h"

using std::map;
using std::string;
using std::list;
using std::pair;
using std::function;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;
using ogdf::edge;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::ArgException;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;

#include <mco/ep/dual_benson/ep_dual_benson.h>
#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/basic/point.h>

using mco::EPDualBensonSolver;
using mco::TemporaryGraphParser;
using mco::Point;

void EpBensonModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Dual Benson to find the Pareto-frontier of the efficient path problem.", ' ', "0.1");
        
        ValueArg<double> epsilon_argument("e", "epsilon", "Epsilon to be used in floating point calculations.", false, 1E-8, "epsilon");
        
        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");
        
        cmd.add(epsilon_argument);
        cmd.add(file_name_argument);
        
        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        double epsilon = epsilon_argument.getValue();
        
        Graph graph;
        EdgeArray<Point> costs(graph);
        unsigned dimension;
        node source, target;
        
        TemporaryGraphParser parser;
        
        parser.getGraph(file_name, graph, costs, dimension, source, target);
        
        EPDualBensonSolver<> solver(epsilon);
        
        auto cost_function = [costs] (edge e) { return &costs[e]; };
        
        solver.Solve(graph, cost_function, source, target);
        
        solutions_.insert(solutions_.begin(),
                          solver.solutions().cbegin(),
                          solver.solutions().cend());
        
    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const list<edge>, const Point>>& EpBensonModule::solutions() {
    return solutions_;
}

string EpBensonModule::statistics() {
    string stats("");
    return stats;
}
