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
#include <mco/benchmarks/temporary_graphs_parser.h>

using mco::Point;
using mco::TemporaryGraphParser;

class TemporaryGraphInstanceFixture
: public ::testing::TestWithParam<tuple<string, unsigned, unsigned>> {
protected:
    const string filename_ = get<0>(GetParam());
    const unsigned expected_number_of_nodes_ = get<1>(GetParam());
    const unsigned expected_number_of_edges_ = get<2>(GetParam());
    
};

TEST_P(TemporaryGraphInstanceFixture, NodeEdgeMatch) {
    Graph graph;
    EdgeArray<Point> costs(graph);
    unsigned dimension;
    node source;
    node target;
    
    TemporaryGraphParser parser;
    
    parser.getGraph(filename_,
                    graph,
                    costs,
                    dimension,
                    source,
                    target);
    
    EXPECT_EQ(expected_number_of_nodes_, graph.numberOfNodes());
    EXPECT_EQ(expected_number_of_edges_, graph.numberOfEdges());
}

INSTANTIATE_TEST_CASE_P(TemporaryGraphParsingTest,
                        TemporaryGraphInstanceFixture,
                        Values(
                               make_tuple(string("../../../instances/ep/graph_1000_1000"),
                                          (unsigned) 318,
                                          (unsigned) 1166)
                               ));
