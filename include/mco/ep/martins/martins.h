#pragma once
/*

 * martins.h
 *
 *  Created on: 12.03.2013
 *      Author: fritz
 */

#ifndef MARTINS_B_H_
#define MARTINS_B_H_

#include <mco/basic/abstract_solver.h>

namespace mco {

class EpSolverMartins : public AbstractSolver<std::list<ogdf::edge>> {

public:
	explicit EpSolverMartins(double epsilon = 0)
    :   epsilon_(epsilon),
        do_value_callback_(false),
        value_callback_([] (Point) {return;}),
        do_path_callback_(false),
        path_callback_([] (std::list<ogdf::node>) {return;}) { }
    
    void Solve(ogdf::Graph& graph,
               std::function<const Point*(ogdf::edge)> weights,
               unsigned dimension,
               ogdf::node source,
               ogdf::node target,
               const Point& absolute_bound,
               std::function<double(ogdf::node, unsigned)> heuristic,
               std::list<Point> first_phase_bounds = std::list<Point>(),
               bool directed = true) {
        
        Solve(graph,
              weights,
              dimension,
              source,
              target,
              absolute_bound,
              heuristic,
              std::list<std::pair<ogdf::NodeArray<Point*>,
                                  ogdf::NodeArray<ogdf::edge>>>(),
              first_phase_bounds,
              directed);
    }
    
    void Solve(ogdf::Graph& graph,
               std::function<const Point*(ogdf::edge)> weights,
               unsigned dimension,
               ogdf::node source,
               ogdf::node target,
               const Point& absolute_bound,
               std::list<std::pair<ogdf::NodeArray<Point*>,
                                   ogdf::NodeArray<ogdf::edge>>> initial_labels,
               std::function<double(ogdf::node, unsigned)> heuristic,
               bool directed = true) {
        
        Solve(graph,
              weights,
              dimension,
              source,
              target,
              absolute_bound,
              heuristic,
              initial_labels,
              std::list<Point>(),
              directed);
        
    }

    
    void Solve(ogdf::Graph& graph,
               std::function<const Point*(ogdf::edge)> weights,
               unsigned dimension,
               ogdf::node source,
               ogdf::node target,
               bool directed = true) {
        
        Point absolute_bound(numeric_limits<double>::infinity(), dimension);
        
        Solve(graph,
              weights,
              dimension,
              source,
              target,
              absolute_bound,
              [] (ogdf::node, unsigned) { return 0; },
              std::list<Point>(),
              directed);
                           
    }
    
    void set_value_callback(std::function<void(Point)> callback) {
        value_callback_ = callback;
        do_value_callback_ = true;
    }
    
    void set_path_callback(std::function<void(std::list<ogdf::node>)> callback) {
        path_callback_ = callback;
        do_path_callback_ = true;
    }
    
private:
    const double epsilon_;
    
    bool do_value_callback_;
    std::function<void(Point)> value_callback_;
    bool do_path_callback_;
    std::function<void(std::list<ogdf::node>)> path_callback_;
    
    void Solve(ogdf::Graph& graph,
               std::function<const Point*(ogdf::edge)> weights,
               unsigned dimension,
               ogdf::node source,
               ogdf::node target,
               const Point& absulute_bound,
               std::function<double(ogdf::node, unsigned)> heuristic,
               std::list<std::pair<ogdf::NodeArray<Point*>,
                                   ogdf::NodeArray<ogdf::edge>>> initial_labels,
               std::list<Point> first_phase_bounds = std::list<Point>(),
               bool directed = true);
    
    struct Label {
        const Point * const point;
        ogdf::node n;
        const Label * const pred;
        bool mark_dominated;
        bool in_queue;
        
        inline Label(const Point *point, ogdf::node n, const Label *pred);
        inline Label(const Label &label);
        
        Label & operator=(const Label &label) = delete;
        
        ~Label() {
            delete point;
        }
    };
    
    void construct_labels(ogdf::NodeArray<std::list<Label*>> & labels,
                          std::list<std::pair<ogdf::NodeArray<Point*>,
                                              ogdf::NodeArray<ogdf::edge>>>& initial_labels,
                          const Point& absolute_bound);
    
    struct LexLabelComp {
        bool operator()(const Label& l1, const Label& l2) {
            return LexPointComparator()(l2.point, l1.point);
        }
        
        bool operator()(const Label* l1, const Label* l2) {
            return LexPointComparator()(l2->point, l1->point);
        }
    };
    
    struct HeuristicLexLabelComp {
        HeuristicLexLabelComp(unsigned dimension,
                              std::function<double(ogdf::node, unsigned)> heuristic)
        :   heuristic_(heuristic),
            dimension_(dimension) {}
        
        bool operator()(const Label& l1, const Label& l2) {
            unsigned i = 0;
            for(; i < dimension_; ++i) {
                if(l1.point->operator[](i) + heuristic_(l1.n, i) <
                   l2.point->operator[](i) + heuristic_(l2.n, i)) {
                    return false;
                } else if(l1.point->operator[](i) + heuristic_(l1.n, i) >
                          l2.point->operator[](i) + heuristic_(l2.n, i)) {
                    return true;
                }
            }
            return false;
        }
        
        bool operator()(const Label* l1, const Label* l2) {
            return operator()(*l1, *l2);
        }
        
    private:
        std::function<double(ogdf::node, unsigned)> heuristic_;
        unsigned dimension_;
    };

};
    
EpSolverMartins::Label::
Label(const Point *point,
      ogdf::node n,
      const Label *pred)
:   point(point),
    n(n),
    pred(pred),
    mark_dominated(false),
    in_queue(true) {
}

EpSolverMartins::Label::
Label(const Label &label)
:   point(new Point(*label.point)),
    n(label.n),
    pred(label.pred),
    mark_dominated(label.mark_dominated),
    in_queue(label.in_queue) {
}

}

#endif /* MARTINS_H_ */
