//
//  main.cpp
//  cli
//
//  Created by Fritz Bökler on 09.04.14.
//
//

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <chrono>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::list;
using std::pair;
using std::make_pair;
using std::chrono::steady_clock;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::exception;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;
using ogdf::edge;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::SwitchArg;
using TCLAP::ValueArg;

#include <mco/basic/point.h>
#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/ep/basic/dijkstra.h>
#include <mco/ep/brum_shier/ep_solver_bs.h>
#include <mco/ep/brum_shier/ep_weighted_bs.h>
#include <mco/ep/brum_shier/ep_weighted_bs.h>
#include <mco/ep/martins/martins.h>
#include <mco/ep/dual_benson/ep_dual_benson.h>
#include <mco/ep/martins/weighted_martins.h>
#include <mco/generic/benson_dual/ove_cdd.h>
#include <mco/generic/benson_dual/ove_node_lists.h>

using mco::Point;
using mco::TemporaryGraphParser;
using mco::EpSolverBS;
using mco::EpWeightedBS;
using mco::EPDualBensonSolver;
using mco::EpSolverMartins;
using mco::EpWeightedMartins;
using mco::Dijkstra;
using mco::DijkstraModes;
using mco::LexPointComparator;
using mco::ParetoDominationPointComparator;
using mco::ComponentwisePointComparator;
using mco::EqualityPointComparator;

#include "basic/modules.h"
#include "modules/ep_benson_module.h"
#include "modules/ep_martins_module.h"

int main(int argc, char** argv) {
    try {
        ModuleFactory module_factory;
        
        EpBensonModule benson_module;
        EpMartinsModule martins_module;
        
        module_factory.add_module("ep-dual-benson", benson_module);
        module_factory.add_module("ep-martins", martins_module);
        
        list<pair<unsigned, BasicModule*>> modules = module_factory.parse_module_list(argc, argv);
        
        if(modules.size() < 1) {
            cout << "Need to specify a module." << endl;
            return 0;
        }
        
        if(modules.size() > 1) {
            cout << "Only one module allowed yet." << endl;
            return 0;
        }
        
        unsigned argument_position = modules.begin()->first;
        BasicModule* choosen_module = modules.begin()->second;
        
        CmdLine cmd("MCO Library Command Line Tool", ' ',  "0.1");
        
        SwitchArg print_frontier_arg("F", "frontier", "Prints the frontier if the problem was feasible", false);
        
        SwitchArg print_solutions_arg("S", "solutions", "Prints the efficient solutions which have been found if the problem was feasible", false);
        
        SwitchArg force_print_all_arg("f", "force-all", "Prints all points/solutions no matter how many", false);
        
        SwitchArg print_verbose_arg("v", "verbose", "Prints output in human readable form.", false);
        
        SwitchArg print_count_arg("c", "count", "Prints the size of the found Pareto-frontier.", false);
        
        SwitchArg print_timing_arg("t", "timing", "Prints timing information", false);
        
        cmd.add(print_frontier_arg);
        cmd.add(print_solutions_arg);
        cmd.add(force_print_all_arg);
        cmd.add(print_count_arg);
        cmd.add(print_verbose_arg);
        cmd.add(print_timing_arg);
        
        cmd.parse(argument_position, argv);
        
        bool print_frontier = print_frontier_arg.getValue();
        bool print_solutions = print_solutions_arg.getValue();
        bool force_print_all = force_print_all_arg.getValue();
        bool print_count = print_count_arg.getValue();
        bool print_verbose = print_verbose_arg.getValue();
        bool print_timing = print_timing_arg.getValue();
        
        choosen_module->perform(argc - argument_position,
                                argv + argument_position);
        
        
        if(print_frontier || print_solutions || print_count) {
            auto ep_algo_module = dynamic_cast<AlgorithmModule<list<edge>>*>(choosen_module);
            
            auto solutions = ep_algo_module->solutions();
            
            cout << solutions.size() << " points" << endl;
            
            if(print_frontier || print_solutions) {
                
                int count = 0;
                auto solution_it = solutions.cbegin();
                while(solution_it != solutions.cend() && (
                      count < 25 || force_print_all)) {
                    auto solution = *solution_it;
                    
                    if(print_frontier) {
                        auto point_it = solution.second.cbegin();
                        while(point_it != solution.second.cend()) {
                            cout << *point_it++ << ", ";
                        }
                    }
                    if(print_solutions) {
                        for(auto edge : solution.first) {
                            cout << edge->source() << ", ";
                        }
                        cout << (*solution.first.rbegin())->target() << ", ";
                    }
                    cout << endl;
                    solution_it++;
                    count++;
                }
                
            }
        }
        
        if(print_timing) {
            cout << "Timining information" << endl;
        }
        
    } catch (exception& e) {
        cout << e.what() << endl;
    }
    
}

/*int main(int argc, char** argv) {
    if(argc != 3) {
        cout << "Usage: " << argv[0] << "<algorithm> <file>" << endl;
    }
    
    string algorithm(argv[1]);
    string filename(argv[2]);

    Graph graph;
    EdgeArray<Point> costs(graph);
    node source;
    node target;
    unsigned dimension;
    
    mco::TemporaryGraphParser parser;
    
    parser.getGraph(filename,
                    graph,
                    costs,
                    dimension,
                    source,
                    target);
    
    cout << "Read a graph with " << graph.numberOfNodes() <<
    " nodes and " << graph.numberOfEdges() << " edges." << endl;
    cout << "Solving..." << endl;
    
    auto cost_function = [&costs] (edge e) {
        return &costs[e];
    };
    
    if (algorithm.compare("dual-benson") == 0) {
        
        EPDualBensonSolver<> solver;
        
        solver.Solve(graph, cost_function, source, target);
        
        cout << solver.solutions().size() << endl;
    
        for(auto p : solver.solutions()) {
            cout << *p << endl;
        }
        
    } else if(algorithm.compare("label-correcting") == 0) {
    
        EpSolverBS solver;
        
        solver.Solve(graph, cost_function, dimension, source, target, false);
        
        cout << solver.solutions().size() << endl;
        
        vector<const Point*> solutions(solver.solutions().cbegin(),
                                       solver.solutions().cend());
        
        std::sort(solutions.begin(), solutions.end(), mco::LexPointComparator());
        
//        mco::ComponentwisePointComparator comp;
//        mco::EqualityPointComparator eq;
//        bool error = false;
//        
//        auto iter1 = solutions.begin();
//        while(iter1 != solutions.end()) {
//            auto iter2 = iter1 + 1;
//            while(iter2 != solutions.end()) {
//                if(comp(*iter1, *iter2) || comp(*iter2, *iter1) || eq(*iter1, *iter2)) {
//                    cout << *(*iter1) << endl;
//                    cout << *(*iter2) << endl;
//                    error = true;
//                }
//                ++iter2;
//            }
//            ++iter1;
//        }
//        
//        cout << error << endl;
        
        for(auto p : solutions) {
            cout << *p << endl;
        }
    } else if(algorithm.compare("w-label-correcting") == 0) {
        EpWeightedBS solver;
        
        solver.Solve(graph, cost_function, dimension, source, target, false);
        
        cout << solver.solutions().size() << endl;
        
        for(auto p : solver.solutions()) {
            cout << *p << endl;
        }
    } else if(algorithm.compare("martins") == 0) {
        EpSolverMartins solver;
        
        solver.Solve(graph, cost_function, dimension, source, target, false);
        
        cout << solver.solutions().size() << endl;
        
        for(auto p : solver.solutions()) {
            cout << *p << endl;
        }
    } else if(algorithm.compare("w-martins") == 0) {
        EpWeightedMartins solver;
        
        solver.Solve(graph, cost_function, dimension, source, target, false);
        
        cout << solver.solutions().size() << endl;
        
        for(auto p : solver.solutions()) {
            cout << *p << endl;
        }
    } else if(algorithm.compare("pre-martins") == 0) {
        Dijkstra<double> sssp_solver;
        
        vector<NodeArray<double>> distances(dimension, graph);
        NodeArray<edge> predecessor(graph);
        
        cout << "calculating heuristic..." << endl;
        
        for(unsigned i = 0; i < dimension; ++i) {
            auto length = [&costs, i] (edge e) {
                return costs[e][i];
            };
            
            sssp_solver.singleSourceShortestPaths(graph,
                                                  length,
                                                  target,
                                                  predecessor,
                                                  distances[i],
                                                  DijkstraModes::Undirected);
        }
        
        auto heuristic = [&distances] (node n, unsigned objective) {
            return distances[objective][n];
        };
        
        Point absolute_bound(numeric_limits<double>::infinity(), dimension);
        for(unsigned i = 0; i < dimension; ++i) {
            absolute_bound[i] = numeric_limits<double>::infinity();
        }
//        bound[dimension - 1] = 1.4 * distances[dimension - 1][source];
        
        EpSolverMartins solver;
        
        cout << "Running Martins algorithm..." << endl;
        
        solver.Solve(graph,
                     cost_function,
                     dimension,
                     source,
                     target,
                     absolute_bound,
                     heuristic,
                     list<Point>(),
                     false);
        
        cout << "Size of the Pareto-frontier: " << solver.solutions().size() << endl;
        
    } else if(algorithm.compare("fpre-martins") == 0) {
        Dijkstra<double> sssp_solver;
        
        vector<NodeArray<double>> distances(dimension, graph);
        NodeArray<edge> predecessor(graph);
        
        cout << "calculating heuristic..." << endl;
        
        for(unsigned i = 0; i < dimension; ++i) {
            auto length = [&costs, i] (edge e) {
                return costs[e][i];
            };
            
            sssp_solver.singleSourceShortestPaths(graph,
                                                  length,
                                                  target,
                                                  predecessor,
                                                  distances[i],
                                                  DijkstraModes::Undirected);
        }
        
        auto heuristic = [&distances] (node n, unsigned objective) {
            return distances[objective][n];
        };
        
        Point absolute_bound(numeric_limits<double>::infinity(), dimension);
        for(unsigned i = 0; i < dimension - 1; ++i) {
            absolute_bound[i] = numeric_limits<double>::infinity();
        }
        absolute_bound[dimension - 1] = 1.4 * distances[dimension - 1][source];
        
        cout << "Running first phase..." << endl;
        
        list<Point> first_phase_bounds;
        {
        
            EPDualBensonSolver<> weighted_solver;
            
            weighted_solver.Solve(graph, cost_function, source, target);
            
            
            for(auto point : weighted_solver.solutions()) {
                first_phase_bounds.push_back(*point);
            }
            
        }
        
        EpSolverMartins solver;
        
        cout << "Running Martins algorithm..." << endl;
        
        solver.Solve(graph,
                     cost_function,
                     dimension,
                     source,
                     target,
                     absolute_bound,
                     heuristic,
                     first_phase_bounds,
                     false);
        
        cout << "Size of the Pareto-frontier: " << solver.solutions().size() << endl;
        
    } else if(algorithm.compare("flabel-martins") == 0) {
        Dijkstra<double> sssp_solver;
        
        vector<NodeArray<double>> distances(dimension, graph);
        NodeArray<edge> predecessor(graph);
        
        cout << "Calculating heuristic... ";
        steady_clock::time_point start = steady_clock::now();
        
        for(unsigned i = 0; i < dimension; ++i) {
            auto length = [&costs, i] (edge e) {
                return costs[e][i];
            };
            
            sssp_solver.singleSourceShortestPaths(graph,
                                                  length,
                                                  target,
                                                  predecessor,
                                                  distances[i],
                                                  DijkstraModes::Undirected);
        }
        
        auto heuristic = [&distances] (node n, unsigned objective) {
            return distances[objective][n];
        };
        
        steady_clock::time_point heuristic_end = steady_clock::now();
        duration<double> heuristic_computation_span
            = duration_cast<duration<double>>(heuristic_end - start);
        
        cout << "Done. (" << heuristic_computation_span.count() << "s)" << endl;
        
        Point absolute_bound(numeric_limits<double>::infinity(), dimension);
        Point relative_bound(numeric_limits<double>::infinity(), dimension);
        
        for(unsigned i = 0; i < dimension; ++i) {
            absolute_bound[i] = numeric_limits<double>::infinity();
        }
        absolute_bound[dimension - 1] = 1.4 * distances[dimension - 1][source];
        
        absolute_bound[0] = 0.01 * absolute_bound[dimension - 1];    // Bündelung
        absolute_bound[1] = 0.01 * absolute_bound[dimension - 1];    // Flughäfen
        absolute_bound[2] = 0.01 * absolute_bound[dimension - 1];    // Tagebau
        absolute_bound[3] = 0.01 * absolute_bound[dimension - 1];    // VSG
        absolute_bound[4] = 0.01 * absolute_bound[dimension - 1];   // Siedlungen
        
        cout << "Running first phase... ";
        
        list<pair<NodeArray<Point *>, NodeArray<edge>>> solutions;
        
        auto callback = [&solutions, &graph, dimension] (NodeArray<Point *>& distances,
                                     NodeArray<edge>& predecessors) {
            
            NodeArray<Point *> new_distances(graph);
            
            for(auto n : distances.graphOf()->nodes) {
                Point* p = new Point(dimension);
                std::copy(distances[n]->cbegin() + 1, distances[n]->cend(),
                          p->begin());
                new_distances[n] = p;
            }
            
            solutions.push_back(make_pair(new_distances, predecessors));
        };
        
        {
            
            EPDualBensonSolver<> weighted_solver;
            
            weighted_solver.Solve(graph, cost_function, source, target, callback);
            
        }
        
        steady_clock::time_point fp_end = steady_clock::now();
        duration<double> fp_computation_span
            = duration_cast<duration<double>>(fp_end - heuristic_end);
        
        cout << "Done. (" << fp_computation_span.count() << "s)" << endl;

        vector<list<node>> paths;
        vector<Point> values;
        
        auto path_callback = [&paths] (list<node> path) {
            paths.push_back(path);
        };
        
        auto value_callback = [&values] (Point p) {
            values.push_back(p);
        };
        
        EpSolverMartins solver;
        
        solver.set_path_callback(path_callback);
        solver.set_value_callback(value_callback);
        
        cout << "Running Martins algorithm... ";
        std::flush(cout);
        
        solver.Solve(graph,
                     cost_function,
                     dimension,
                     source,
                     target,
                     absolute_bound,
                     solutions,
                     heuristic,
                     false);
        
        steady_clock::time_point martins_end = steady_clock::now();
        duration<double> martins_computation_span
            = duration_cast<duration<double>>(martins_end - fp_end);
        
        cout << "Done. (" << martins_computation_span.count() << "s)" << endl;

        cout << "Size of the Pareto-frontier: " << solver.solutions().size() << endl;
        
        if(values.size() < 25) {
            for(auto value : values) {
                cout << value << endl;
            }
        }
        
    } else {
        cout << "Unknown algorithm: " << algorithm << endl;
    }
    
    
}*/