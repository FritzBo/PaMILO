//
//  ap_hungarian_benson_test.cpp
//  mco
//
//  Created by Fritz Bökler on 06.04.14.
//
//

#include <set>
#include <tuple>
#include <string>

using std::set;
using std::tuple;
using std::string;
using std::get;
using std::make_tuple;

#include <gtest/gtest.h>

using ::testing::Values;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::node;
using ogdf::EdgeArray;

#include <mco/basic/point.h>
#include <mco/ap/basic/ap_instance.h>
#include <mco/ap/benson_dual/ap_benson_dual_solver.h>
#include <mco/benchmarks/mcap_parser.h>

using mco::Point;
using mco::MCAPParser;
using mco::AssignmentInstance;
using mco::APBensonDualSolver;
using mco::EqualityPointComparator;

class ParetoInstanceTestFixture
: public ::testing::TestWithParam<tuple<string, unsigned>> {
protected:
    const string filename_ = get<0>(GetParam());
    const unsigned expected_number_of_minimizers_ = get<1>(GetParam());
    
};

TEST_P(ParetoInstanceTestFixture, CountMatch) {
    Graph graph;
    EdgeArray<Point*> costs(graph);
    set<node> agents;
    
    MCAPParser parser(filename_);
    
    AssignmentInstance instance = parser.get_instance(graph,
                                                      costs,
                                                      agents);
    
    APBensonDualSolver<> solver;
    
    solver.Solve(instance);
    
    EXPECT_EQ(expected_number_of_minimizers_, solver.solutions().size());
        
    for(auto e : graph.edges) {
        delete costs(e);
    }
}

INSTANTIATE_TEST_CASE_P(InstanceTests,
                        ParetoInstanceTestFixture,
                        Values(
                               make_tuple(string("../../../instances/ap/pge_3_10_2"), 28),
                               make_tuple(string("../../../instances/ap/pge_3_05_5"), 12),
                               make_tuple(string("../../../instances/ap/map_1_2_0"), 1)
                               ));
