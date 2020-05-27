/**
 * dijkstra.cc
 *
 *  Created on: 22.03.2014
 *      Author: Denis Kurz
 */

#include <functional>
#include <limits>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>
#include <mco/basic/lex_point_comparator.h>
#include <mco/ep/basic/binary_heap.h>
#include <mco/ep/basic/dijkstra.h>

using std::function;
using std::numeric_limits;

using ogdf::Graph;
using ogdf::node;
using ogdf::edge;
using ogdf::NodeArray;
using ogdf::EdgeArray;

namespace mco {

void LexDijkstra::singleSourceShortestPaths(
        Graph const &graph,
        function<Point*(edge)> weight,
        node const source,
        NodeArray<Point *> &distance,
        NodeArray<edge> &predecessor,
        function<bool(node,edge)> mode) {
    
    LexPointComparator less;
    
    unsigned dim = weight(graph.chooseEdge())->dimension();

    BinaryHeap2<Point *, node> queue(graph.numberOfNodes(), less);
    NodeArray<int> qpos(graph);
    for(auto v : graph.nodes) {
        for(unsigned i = 0; i < dim; ++i) {
            (*distance[v])[i] = numeric_limits<double>::max();
        }
        predecessor[v] = nullptr;
        queue.insert(v, distance[v], &qpos[v]);
    }

    *distance[source] = Point(0.0, distance[source]->dimension());
    
    queue.decreaseKey(qpos[source], distance[source]);

    Point tmp(dim);
    while (!queue.empty()) {
        node v = queue.extractMin();
        for(auto adj : v->adjEdges) {
            edge e = adj->theEdge();
            if(!mode(v, e)) continue;
            node w = e->opposite(v);
            tmp = *distance[v];
            tmp += *weight(e);
            if(less(&tmp, distance[w])) {
                *distance[w] = tmp;
                queue.decreaseKey(qpos[w], distance[w]);
                predecessor[w] = e;
            }
        }
    }
}

function<bool(node,edge)> const DijkstraModes::Forward =
    [](node v, edge e) {
        return v == e->source();
    };

function<bool(node,edge)> const DijkstraModes::Backward =
    [](node v, edge e) {
        return v == e->target();
    };

function<bool(node,edge)> const DijkstraModes::Undirected =
    [](node, edge) {
        return true;
    };

}

