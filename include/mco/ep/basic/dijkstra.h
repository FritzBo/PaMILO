/**
 * dijkstra.h
 *
 *  Created on: 25.11.2013
 *      Author: Denis Kurz
 */

#ifndef MCO_DIJKSTRA_H_
#define MCO_DIJKSTRA_H_

#include <functional>
#include <limits>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>
#include <mco/ep/basic/binary_heap.h>

namespace mco {

//! TODO doxygen
struct DijkstraModes {

    // TODO doxygen
    static std::function<bool(ogdf::node,ogdf::edge)> const Forward;

    // TODO doxygen
    static std::function<bool(ogdf::node,ogdf::edge)> const Backward;

    // TODO doxygen
    static std::function<bool(ogdf::node,ogdf::edge)> const Undirected;

};

//! TODO doxygen
template<typename T>
class Dijkstra {

public:

    //! TODO doxygen
    void singleSourceShortestPaths(
            ogdf::Graph const &graph,
            std::function<T(ogdf::edge)> weight,
            ogdf::node const source,
            ogdf::NodeArray<ogdf::edge> &predecessor,
            ogdf::NodeArray<T> &distance,
            std::function<bool(ogdf::node,ogdf::edge)> mode=DijkstraModes::Forward) {
        // Initialization: set distances
        for(auto v : graph.nodes) {
            distance[v] = numeric_limits<T>::max();
        }

        // Initialization: populate priority queue
        mco::BinaryHeap2<T, ogdf::node> queue(graph.numberOfNodes());
        ogdf::NodeArray<int> qpos(graph);
        for(auto v : graph.nodes) {
            queue.insert(v, distance[v], &qpos[v]);
        }
        distance[source] = 0;
        queue.decreaseKey(qpos[source], distance[source]);

        // Dijkstra: empty queue, update distances accordingly
        while(!queue.empty()) {
            auto v = queue.extractMin();
            for(auto adj : v->adjEdges) {
                auto e = adj->theEdge();
                if(!mode(v, e)) continue;
                auto w = e->opposite(v);
                T newDist = distance[v] + weight(e);
                if(distance[w] > newDist) {
                    queue.decreaseKey(qpos[w], (distance[w] = newDist));
                    predecessor[w] = e;
                }
            }
        }
    }

};


//! TODO doxygen
class LexDijkstra {

public:

    //! TODO doxygen
    void singleSourceShortestPaths(
            ogdf::Graph const &graph,
            std::function<Point*(ogdf::edge)> weight,
            ogdf::node const source,
            ogdf::NodeArray<Point *> &distance,
            ogdf::NodeArray<ogdf::edge> &predecessor,
            std::function<bool(ogdf::node,ogdf::edge)> mode=DijkstraModes::Forward);

};

}

#endif /* MCO_DIJKSTRA_H_ */
