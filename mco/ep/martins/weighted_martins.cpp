/*
 * martins.cpp
 *
 *  Created on: 23.03.2013
 *      Author: fritz
 */

#include <mco/ep/martins/weighted_martins.h>

#include <queue>
#include <vector>
#include <set>
#include <list>
#include <iomanip>

using std::priority_queue;
using std::vector;
using std::set;
using std::list;
using std::function;
using std::pair;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::AdjElement;
using ogdf::node;
using ogdf::edge;
using ogdf::NodeArray;

#include <mco/basic/point.h>
#include <mco/ep/basic/ep_instance.h>
#include <mco/ep/martins/label.h>

namespace mco {
    
bool EpWeightedMartins::NodeEntry::
add_label(Label& new_label) {
    
//    cout << "new label:" << endl;
//    cout << *new_label.point << endl;
//    cout << "old labels:" << endl;
//    for(auto label : label_set) {
//        if(label != nullptr) {
//            cout << *label->point << endl;
//        }
//    }
//    cout << endl;
    
    // If hull matrix is not set, set it
    if(hull_matrix == nullptr) {
        initialHullMatrix(dimension_, new_label);
        return true;
    
    // If there are only two points to compare
    } else if(head ==  dimension_ + 1) {
        ComponentwisePointComparator cmp(0, false);
        auto point1 = *label_set[0]->point;
        auto point2 = *new_label.point;
        
        if(cmp(point1, point2)) {
            return false;
        } else if(cmp(point2, point1)) {
            label_set[0] = &new_label;
            
            for(unsigned i = 0; i < dimension_; ++i) {
                dd_set_d(hull_matrix->matrix[head - 1][i + 1], point2[i]);
            }
            dd_set_d(hull_matrix->matrix[head - 1][0], 1.0);
            
            return true;

        } else {
            
            for(unsigned i = 0; i < dimension_; ++i) {
                dd_set_d(hull_matrix->matrix[head][i + 1], point2[i]);
            }
            dd_set_d(hull_matrix->matrix[head][0], 1.0);
            
            label_set[1] = &new_label;
            
            ++head;
            
            return true;
        }
    } else {
        
        for(auto label : label_set) {
            if(label != nullptr) {
                if(ComponentwisePointComparator()(label->point, new_label.point)) {
                    return false;
                }
            }
        }
    
        // If we have not enough space, extend hull matrix
        if(hull_matrix->rowsize <= head && free_list.empty()) {
            extend_hull_matrix();
        }
        
        const Point& new_cost = *new_label.point;
        
        // Find position to insert and insert point in hull_matrix at
        // a free position or at head pointer
        dd_rowrange position;
        bool used_free_list = false;
        if(free_list.empty()) {
            position = head;
        } else {
            position = free_list.front();
            used_free_list = true;
        }
        
        for(unsigned i = 0; i < dimension_; ++i) {
            dd_set_d(hull_matrix->matrix[position][i + 1], new_cost[i]);
        }
        dd_set_d(hull_matrix->matrix[position][0], 1.0);
        
//        dd_WriteMatrix(stdout, hull_matrix);
        
        // Check if new point is redundant
        dd_Arow cert = new double [dimension_ + 1][1];
        dd_ErrorType err = dd_NoError;
        if(dd_Redundant(hull_matrix, position + 1, cert, &err)) {
            
            for(unsigned i = 0; i < dimension_ + 1; ++i) {
                dd_set_d(hull_matrix->matrix[position][i], 0.0);
            }
            
            return false;
        }
        
        label_set[position - dimension_] = &new_label;
        
        if(used_free_list) {
            free_list.pop_front();
        } else {
            ++head;
        }
        
        // Check if other points are now redundant
        for(unsigned i = dimension_; i < head; ++i) {
            Label* label = label_set[i - dimension_];
            if(label != nullptr && label->in_queue && i != position) {
                if(dd_Redundant(hull_matrix, i + 1, cert, &err)) {
                    
                    for(unsigned j = 0; j < dimension_ + 1; ++j) {
                        dd_set_d(hull_matrix->matrix[i][j], 0.0);
                    }
                    
                    if(i < head - 1) {
                        free_list.push_back(i);
                    } else {
                        --head;
                    }
                 
                    label->mark_dominated = true;
                    label_set[i - dimension_] = nullptr;
                }
            }
        }
        
//        cout << "new label set:" << endl;
//        for(auto label : label_set) {
//            if(label != nullptr) {
//                cout << *label->point << endl;
//            }
//        }
//        cout << endl;
//        cout << endl;
        
        return true;
    }
}
    
void EpWeightedMartins::
Solve(Graph& graph,
      function<const Point*(edge)> weights,
      unsigned dimension,
      node source,
      node target,
      bool directed) {
    
    using LabelPriorityQueue = priority_queue<Label *, vector<Label *>, LexLabelComp>;
    
    dd_set_global_constants();
    
	LabelPriorityQueue lex_min_label((LexLabelComp()));
	NodeArray<NodeEntry> node_entry(graph, dimension);

	Label *null_label = new Label(Point::Null(dimension), source, nullptr);
    null_label->in_queue = true;
    node_entry[source].label_set.push_back(null_label);
    
	lex_min_label.push(null_label);

	while(!lex_min_label.empty()) {
		Label *label = lex_min_label.top();
		lex_min_label.pop();
        assert(label->in_queue);
        label->in_queue = false;

		if(label->mark_dominated) {
			delete label;
			continue;
		}

		const Point *label_cost = label->point;
		node n = label->n;
        
        if(n == target) {
//            cout << *label_cost << endl;
            continue;
        }

//		cout << endl << n << ", " << *label->point << ": ";

		AdjElement *adj;
		forall_adj(adj, n) {
			edge e = adj->theEdge();

			if(e->isSelfLoop())
				continue;

			node v = e->target();

            if(directed) {
                if(v == n)
                    continue;
            } else {
                if(v == n) {
                    v = e->source();
                }
            }

			if(v == source)
				continue;

//			cout << v << ", ";

			const Point *edge_cost = weights(e);
			const Point *new_cost = new Point(*edge_cost + *label_cost);			// Owner

			Label * new_label = new Label(new_cost, v, label);
            
            if(node_entry[v].add_label(*new_label)) {
                lex_min_label.push(new_label);
                new_label->in_queue = true;
            } else {
                delete new_label;
            }
        }

	}
    
	list<pair<const list<edge>, const Point>> solutions;
    
	for(auto label : node_entry[target].label_set)
        if(label != nullptr) {
            list<edge> path;
            const Label* curr = label;
            while(curr->n != source) {
                for(auto adj: curr->n->adjEdges) {
                    edge e = adj->theEdge();
                    if(e->source() == curr->pred->n && e->target() == curr->n) {
                        path.push_back(e);
                        break;
                    }
                }
                curr = curr->pred;
            }
            
            path.reverse();
            
            solutions.push_back(make_pair(path, *label->point));
        }

	reset_solutions();
    
	add_solutions(solutions.begin(), solutions.end());

//	for(Label *label : labels[instance().target()]) {
//		const Label *current_label = label;
//		cout << "(";
//		while(current_label->n != instance().source()) {
//			cout << current_label->n << ", ";
//			current_label = current_label->pred;
//		}
//		cout << current_label->n << ")" << endl;
//	}

	node n;
	forall_nodes(n, graph) {
		for(auto &label : node_entry[n].label_set)
			delete label;
	}
    
    dd_free_global_constants();
}

}
