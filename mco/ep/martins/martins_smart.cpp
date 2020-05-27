/*
 * martins.cpp
 *
 *  Created on: 12.03.2013
 *      Author: fritz
 */

#include <mco/ep/martins/martins_smart.h>

#include <queue>
#include <memory>
#include <vector>
#include <set>
#include <list>

using std::priority_queue;
using std::shared_ptr;
using std::make_shared;
using std::vector;
using std::set;
using std::pair;
using std::list;

#include <ogdf/basic/Graph.h>

using ogdf::AdjElement;
using ogdf::node;
using ogdf::edge;
using ogdf::NodeArray;

#include <mco/ep/basic/ep_instance.h>
#include <mco/basic/point.h>

namespace mco {

EpSolverMartinsSmart::EpSolverMartinsSmart(EpInstance &instance): AbstractEpSolver(instance) {

}

struct LabelSmart {
	shared_ptr<const Point> point;
	node n;
	shared_ptr<const LabelSmart> pred;
	bool mark_dominated;

	LabelSmart(shared_ptr<const Point> point, node n, shared_ptr<const LabelSmart> pred) : point(point), n(n), pred(pred), mark_dominated(false) {}
};

// TODO: Functor or lambda
bool lexicographic_smaller_label_smart(shared_ptr<const LabelSmart> p1, shared_ptr<const LabelSmart> p2) {
	return LexPointComparator::is_lex_le(*p1->point, *p2->point, 0);
}

using LabelPriorityQueue    = priority_queue<shared_ptr<const LabelSmart>,
                                             vector<shared_ptr<const LabelSmart>>,
                                             decltype(&lexicographic_smaller_label_smart)>;
    
using LabelSet              = set<shared_ptr<const LabelSmart>, decltype(&lexicographic_smaller_label_smart)>;
using LabelList             = list<shared_ptr<const LabelSmart>>;
using LabelListNodeArray    = NodeArray<list<shared_ptr<LabelSmart>>>;

void EpSolverMartinsSmart::Solve() {
	LabelPriorityQueue lex_min_label(&lexicographic_smaller_label_smart);
	LabelListNodeArray labels(instance().graph());

	shared_ptr<const LabelSmart> null_label = make_shared<LabelSmart>(shared_ptr<const Point>(Point::Null(instance().dimension())), instance().source(), nullptr);
	lex_min_label.push(null_label);
    
    ComponentwisePointComparator comp_leq(0.0, false);

	while(!lex_min_label.empty()) {
		shared_ptr<const LabelSmart> label = lex_min_label.top();
		lex_min_label.pop();

		if(label->mark_dominated)
			continue;

		const shared_ptr<const Point> &label_cost = label->point;
		node n = label->n;

		AdjElement *adj;
		forall_adj(adj, n) {
			edge e = adj->theEdge();

			if(e->isSelfLoop())
				continue;

			node v = e->target();

			if(v == n)
				continue;

			if(v == instance().source())
				continue;

			auto edge_cost = instance().weights()(e);
			auto new_cost = make_shared<const Point>(*edge_cost + *label_cost);

			list<shared_ptr<LabelSmart> > nondominatedLabels;
			bool dominated = false;

			for(auto &target_label : labels[v]) {

				if(comp_leq(*target_label->point, *new_cost)) {
					dominated = true;
					break;
				}

				if(comp_leq(*new_cost, *target_label->point)) {
					target_label->mark_dominated = true;
				} else {
					nondominatedLabels.emplace_back(target_label);
                }
			}

			if(dominated)
				continue;

			shared_ptr<LabelSmart> new_label = make_shared<LabelSmart>(new_cost, v, label);
			labels[v] = nondominatedLabels;
			labels[v].emplace_back(new_label);

			if(v != instance().target())
				lex_min_label.push(new_label);
		}

	}

//	for(auto label : labels[instance().target()]) {
//		auto current_label = label;
//		cout << "(";
//		while(current_label->n != instance().source()) {
//			cout << current_label->n << ", ";
//			current_label = current_label->pred;
//		}
//		cout << current_label->n << ")" << endl;
//	}
    
    list<pair<const list<edge>, const Point>> solutions;
    
	for(auto label : labels[instance().target()])
        if(label != nullptr) {
            list<edge> path;
            shared_ptr<const LabelSmart> curr = label;
            while(curr->n != instance().source()) {
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
}

}
