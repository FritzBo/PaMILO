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

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::SwitchArg;
using TCLAP::ValueArg;

#include <mco/basic/point.h>
#include <mco/pilp/pilp_dual_benson.h>
#include <mco/generic/benson_dual/ove_cdd.h>
#include <mco/generic/benson_dual/ove_node_lists.h>

using mco::Point;
//using mco::EpSolverBS;
//using mco::EpWeightedBS;
//using mco::EpSolverMartins;
//using mco::EpWeightedMartins;
using mco::LexPointComparator;
using mco::ParetoDominationPointComparator;
using mco::ComponentwisePointComparator;
using mco::EqualityPointComparator;

#include "basic/modules.h"
#include "modules/pilp_benson_module.h"

int main(int argc, char** argv) {
try {
	ModuleFactory module_factory;
	
	PilpBensonModule pilp_benson_module;
	
	module_factory.add_module("pilp-dual-benson", pilp_benson_module);
	
	list<pair<unsigned, BasicModule*>> modules = module_factory.parse_module_list(argc, argv);
	
	if(modules.size() < 1) {
		std::cout << "Need to specify a module." << std::endl;
		return 0;
	}
	
	if(modules.size() > 1) {
		std::cout << "Only one module allowed yet." << std::endl;
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
	
	
	/*
	if(print_frontier || print_solutions || print_count) {
		auto ep_algo_module = dynamic_cast<AlgorithmModule<list<edge>>*>(choosen_module);
		
		auto solutions = ep_algo_module->solutions();
		
		std::cout << solutions.size() << " points" << std::endl;
		
		if(print_frontier || print_solutions) {
			
			int count = 0;
			auto solution_it = solutions.cbegin();
			while(solution_it != solutions.cend() && (
				  count < 25 || force_print_all)) {
				auto solution = *solution_it;
				
				if(print_frontier) {
					auto point_it = solution.second.cbegin();
					while(point_it != solution.second.cend()) {
						std::cout << *point_it++ << ", ";
					}
				}
				if(print_solutions) {
					for(auto edge : solution.first) {
						std::cout << edge->source() << ", ";
					}
					std::cout << (*solution.first.rbegin())->target() << ", ";
				}
				std::cout << std::endl;
				solution_it++;
				count++;
			}
			
		}
	}
	*/
	
	if(print_timing) {
		std::cout << "Timining information" << std::endl;
	}
	
} catch (exception& e) {
	std::cout << e.what() << std::endl;
}

}

