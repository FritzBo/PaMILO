//
//  temporary_graphs_parser.h
//  mco
//
//  Created by Fritz Bökler on 08.04.14.
//
//

#ifndef __mco__temporary_graphs_parser__
#define __mco__temporary_graphs_parser__

#include <string>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {
    
class TemporaryGraphParser {
public:
    void getGraph(string filename,
                  ogdf::Graph &graph,
                  ogdf::EdgeArray<mco::Point> &weights,
                  unsigned& dimension,
                  ogdf::node& source,
                  ogdf::node& target,
                  bool cut_precision = false,
                  unsigned precision = 7);
    
private:
    unsigned precision_;
};
    
}

#endif /* defined(__mco__temporary_graphs_parser__) */
