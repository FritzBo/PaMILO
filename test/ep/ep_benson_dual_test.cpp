//
//  ep_benson_dual_test.cpp
//  mco
//
//  Created by Fritz Bökler on 10.04.14.
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
using ogdf::edge;
using ogdf::EdgeArray;

#include <mco/basic/point.h>
#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/ep/dual_benson/ep_dual_benson.h>
#include <mco/generic/benson_dual/ove_cdd.h>

using mco::Point;
using mco::EPDualBensonSolver;
using mco::TemporaryGraphParser;
using mco::OnlineVertexEnumeratorCDD;


class ParetoInstanceTestFixture
: public ::testing::TestWithParam<tuple<string, unsigned>> {
protected:
    const string filename_ = get<0>(GetParam());
    const unsigned expected_number_of_minimizers_ = get<1>(GetParam());
    
};

TEST_P(ParetoInstanceTestFixture, CountMatch) {
    Graph graph;
    EdgeArray<Point> costs(graph);
    unsigned dimension;
    node source;
    node target;
    
    TemporaryGraphParser parser;
    
    parser.getGraph(filename_, graph, costs, dimension, source, target);
    
    EPDualBensonSolver<OnlineVertexEnumeratorCDD> solver;
    
    auto weight_function = [costs] (edge e) {
        return &costs(e);
    };
    
    solver.Solve(graph,
                 weight_function,
                 source,
                 target);

#ifndef NDEBUG
    std::cout << solver.solutions().size() << std::endl;
    
    for(auto p : solver.solutions()) {
        std::cout << p.second << std::endl;
    }
#endif
    
    EXPECT_EQ(expected_number_of_minimizers_, solver.solutions().size());
}

INSTANTIATE_TEST_CASE_P(InstanceTests,
                        ParetoInstanceTestFixture,
                        Values(
//                               make_tuple(string("../../../instances/ep/graph_1000_1000"), (unsigned) 1421),
                               make_tuple(string("../../../instances/ep/grid50_1_1"), (unsigned) 2),
                               make_tuple(string("../../../instances/ep/grid50_50_7"), (unsigned) 13)
                               ));
