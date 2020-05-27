#pragma once
/*

 * martins.h
 *
 *  Created on: 12.03.2013
 *      Author: fritz
 */

#ifndef WEIGHTED_MARTINS_B_H_
#define WEIGHTED_MARTINS_B_H_

#include <mco/basic/abstract_solver.h>

#include <setoper.h>
#include <cdd.h>

namespace mco {

class EpWeightedMartins : public AbstractSolver<std::list<ogdf::edge>> {

public:
	explicit EpWeightedMartins(double epsilon = 0)
    :   comp_leq_(epsilon, false) { }
    
	virtual void Solve(ogdf::Graph& graph,
                       std::function<const Point*(ogdf::edge)> weights,
                       unsigned dimension,
                       ogdf::node source,
                       ogdf::node target,
                       bool directed = true);
    
private:
    ComponentwisePointComparator comp_leq_;
    
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
    
    struct NodeEntry {
        std::vector<Label*> label_set;
        std::list<unsigned> free_list;
        unsigned head;
        dd_MatrixPtr hull_matrix;
        unsigned dimension_;
        
        inline NodeEntry(unsigned dimension)
        :   hull_matrix(nullptr),
            dimension_(dimension) { }
        
        inline bool add_label(Label& new_label);
        
    private:
        inline void initialHullMatrix(unsigned dimension,
                                      Label& new_label);
        inline void extend_hull_matrix();
    };
    
    struct LexLabelComp {
        bool operator()(const Label& l1, const Label& l2) {
            return LexPointComparator()(l2.point, l1.point);
        }
        
        bool operator()(const Label* l1, const Label* l2) {
            return LexPointComparator()(l2->point, l1->point);
        }
    };

};
    
EpWeightedMartins::Label::
Label(const Point *point,
      ogdf::node n,
      const Label *pred)
:   point(point),
    n(n),
    pred(pred),
    mark_dominated(false),
    in_queue(false) {
}

EpWeightedMartins::Label::
Label(const Label &label)
:   point(new Point(*label.point)),
    n(label.n),
    pred(label.pred),
    mark_dominated(label.mark_dominated),
    in_queue(label.in_queue) {
}
    
void EpWeightedMartins::NodeEntry::
initialHullMatrix(unsigned dimension, Label& new_label) {
    hull_matrix = dd_CreateMatrix(dimension * 2,
                                  dimension + 1);
    
    for(unsigned i = 0; i < dimension; ++i) {
        for(unsigned j = 0; j < dimension + 1; ++j) {
            dd_set_d(hull_matrix->matrix[i][j], i + 1 == j ? 1.0 : 0.0);
        }
        dd_set_d(hull_matrix->matrix[i][0], 0);
    }
    
    for(unsigned j = 0; j < dimension; ++j) {
        dd_set_d(hull_matrix->matrix[dimension][j + 1], new_label.point->operator[](j));
    }
    dd_set_d(hull_matrix->matrix[dimension][0], 1.0);
    
    hull_matrix->representation = dd_Generator;
    
    head = dimension + 1;
    
    label_set.resize(dimension);
    label_set[0] = &new_label;
    
}
    
void EpWeightedMartins::NodeEntry::
extend_hull_matrix() {
    unsigned new_rowsize = 2 * hull_matrix->rowsize;
    unsigned colsize = hull_matrix->colsize;
    dd_MatrixPtr new_matrix = dd_CreateMatrix(new_rowsize, colsize);
    
    for(unsigned i = 0; i < hull_matrix->rowsize; ++i) {
        for(unsigned j = 0; j < hull_matrix->colsize; ++j) {
            double old_value = dd_get_d(hull_matrix->matrix[i][j]);
            dd_set_d(new_matrix->matrix[i][j], old_value);
        }
    }
    
    new_matrix->representation = dd_Generator;
    
    dd_FreeMatrix(hull_matrix);
    
    hull_matrix = new_matrix;
    
    label_set.resize(new_rowsize - dimension_);
}



}   // namespace mco

#endif /* WEIGHTED_MARTINS_H_ */
