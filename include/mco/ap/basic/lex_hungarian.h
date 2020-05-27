//
//  hungarian.h
//  mco
//
//  Created by Fritz Bökler on 03.04.14.
//
//

#ifndef __mco__lex__hungarian__
#define __mco__lex__hungarian__

#include <functional>
#include <list>
#include <deque>
#include <set>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

class LexHungarianMethod {
public:
    
    template<typename ConstIterator>
    inline Point Solve(const ogdf::Graph& graph,
                       std::function<Point* (ogdf::edge)> edge_costs,
                       unsigned dimension,
                       ConstIterator agents_begin,
                       ConstIterator agents_end);
    
    Point Solve(const ogdf::Graph& graph,
                std::function<Point* (ogdf::edge)> edge_costs,
                unsigned dimension,
                const std::set<ogdf::node>& agents);

    
    Point value() {
        return value_;
    }
    
private:
    Point value_;
    
};
    
template<typename ConstIterator>
inline Point LexHungarianMethod::
Solve(const ogdf::Graph& graph,
      std::function<Point* (ogdf::edge)> edge_costs,
      unsigned dimension,
      ConstIterator agents_begin,
      ConstIterator agents_end) {

    std::set<ogdf::node> agents(agents_begin, agents_end);
    Solve(graph, edge_costs, dimension, agents);
}

    
}

#endif /* defined(__mco__hungarian__) */
