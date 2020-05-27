//
//  lex_hungarian_test.cpp
//  mco
//
//  Created by Fritz Bökler on 04.04.14.
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
#include <mco/ap/basic/lex_hungarian.h>
#include <mco/benchmarks/mcap_parser.h>

using mco::Point;
using mco::MCAPParser;
using mco::AssignmentInstance;
using mco::LexHungarianMethod;
using mco::EqualityPointComparator;

class InstanceTestFixture
: public ::testing::TestWithParam<tuple<string, Point>> {
protected:
    const string filename_ = get<0>(GetParam());
    const Point expected_outcome_ = get<1>(GetParam());
    
};

TEST_P(InstanceTestFixture, PointMatch) {
    Graph graph;
    EdgeArray<Point*> costs(graph);
    set<node> agents;
    
    MCAPParser parser(filename_);
    
    parser.get_instance(graph,
                        costs,
                        agents);
    
    LexHungarianMethod solver;
    
    Point p = solver.Solve(graph,
                           costs,
                           costs[graph.firstEdge()]->dimension(),
                           agents);

    EqualityPointComparator eq;
    EXPECT_TRUE(eq(p, expected_outcome_));
    
    std::cout << p << std::endl;
    
    for(auto e : graph.edges) {
        delete costs(e);
    }
}

INSTANTIATE_TEST_CASE_P(InstanceTests,
                        InstanceTestFixture,
                        Values(
                               make_tuple(string("../../../instances/ap/pge_3_50_10"),
                                          Point({9, 458, 478})),
                               make_tuple(string("../../../instances/ap/pge_3_40_2"),
                                          Point({8, 344, 406})),
                               make_tuple(string("../../../instances/ap/pge_3_05_5"),
                                          Point({19, 69, 63})),
                               make_tuple(string("../../../instances/ap/map_1_2_0"),
                                          Point({0}))
                        ));
