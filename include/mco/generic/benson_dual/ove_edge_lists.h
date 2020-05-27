#pragma once
/*
 * online_vertex_enumerator.h
 *
 *  Created on: 30.09.2013
 *      Author: fritz
 */

#ifndef EDGE_LIST_VE_H_
#define EDGE_LIST_VE_H_

#include <vector>
#include <map>
#include <queue>
#include <list>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>
#include <mco/generic/benson_dual/abstract_online_vertex_enumerator.h>

namespace mco {

class EdgeListVE : protected AbstractOnlineVertexEnumerator {
public:
	EdgeListVE() = delete;
	virtual inline ~EdgeListVE();

	bool has_next();
	Point * next_vertex();
	void add_hyperplane(Point &vertex, Point &normal, double rhs);

	unsigned int number_of_hyperplanes() {
        return list_of_inequalities_.size();
    }

	EdgeListVE(Point &initial_value, unsigned int dimension, double epsilon);


protected:
	ogdf::Graph vertex_graph_;
	ogdf::NodeArray<Point *> node_points_;

	LexPointComparator comp_;
    EqualityPointComparator eq_comp_;
    
	std::map<Point *, ogdf::node, LexPointComparator> point_nodes_;
	std::priority_queue<Point *, std::vector<Point *>, LexPointComparator> unprocessed_projective_points_;

	std::vector<Point *> list_of_inequalities_;
	ogdf::NodeArray<std::list<int> *> node_inequality_indices_;
	ogdf::NodeArray<unsigned int> birth_index_;
	ogdf::EdgeArray<std::list<int> *> edge_inequality_indices_;

	ogdf::node get_node(Point &non_projective_point);
	Point * to_projective(Point &non_projective_point);
	Point normalize_projective(Point projective_point);

	bool inside_face(ogdf::node n1, ogdf::node n2, bool nondegenerate);
};
    
inline EdgeListVE::
~EdgeListVE() {
    for(auto n : vertex_graph_.nodes) {
        delete node_points_[n];
    }
    
    for(auto inequality : list_of_inequalities_) {
        delete inequality;
    }
}


} /* namespace mco */
#endif /* EDGE_LIST_VE_H_ */
