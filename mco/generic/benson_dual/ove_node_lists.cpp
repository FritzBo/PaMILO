/*
 * online_vertex_enumerator.cpp
 *
 *  Created on: 30.09.2013
 *      Author: fritz
 */

#include <mco/generic/benson_dual/ove_node_lists.h>

#include <cassert>
#include <set>
#include <list>
#include <vector>
#include <cmath>

using std::make_pair;
using std::vector;
using std::list;
using std::set;
using std::abs;

#include <ogdf/basic/Graph.h>

using ogdf::node;
using ogdf::edge;
using ogdf::NodeArray;
using ogdf::AdjElement;

namespace mco {

NodeListVE::
NodeListVE(const Point& initial_value,
           unsigned dimension,
           double epsilon)
:   AbstractOnlineVertexEnumerator(dimension, epsilon),
    node_points_(vertex_graph_, nullptr),
    comp_(epsilon),
    point_nodes_(comp_),
    unprocessed_projective_points_(LexPointComparator(epsilon, false)),
    node_inequality_indices_(vertex_graph_, nullptr),
    birth_index_(vertex_graph_, -1) {
    
    node n, v;
    Point *p;
    for(unsigned int i = 0; i < dimension_ - 1; ++i) {
        n = vertex_graph_.newNode();
        node_inequality_indices_[n] = new list<int>();
        p = new Point(dimension_ + 1);
        for(unsigned int j = 0; j < dimension_ - 1; ++j) {
            (*p)[j] = i == j ? 1 : 0;
            
            if(i != j)
                node_inequality_indices_[n]->push_back(j);
        }
        (*p)[dimension_ - 1] = initial_value[i];
        (*p)[dimension_] = 1;
        
        node_inequality_indices_[n]->push_back(dimension_ - 1);
        node_inequality_indices_[n]->push_back(dimension_);
        birth_index_[n] = dimension_;
        
        point_nodes_.insert(make_pair(p, n));
        node_points_[n] = p;
        unprocessed_projective_points_.push(p);
        
        p = new Point(dimension_ + 1);
        for(unsigned int j = 0; j < dimension_; ++j)
            (*p)[j] = i == j ? 1 : 0;
        (*p)[dimension_] = 0;
        
        list_of_inequalities_.push_back(p);
        
        forall_nodes(v, vertex_graph_) {
            if(v != n) {
                vertex_graph_.newEdge(v, n);
                vertex_graph_.newEdge(n, v);
            }
        }
    }
    
    p = new Point(dimension_ + 1);
    for(unsigned int j = 0; j < dimension_ - 1; ++j)
        (*p)[j] = -1;
    (*p)[dimension_ - 1] = 0;
    (*p)[dimension_] = 1;
    
    list_of_inequalities_.push_back(p);
    
    n = vertex_graph_.newNode();
    node_inequality_indices_[n] = new list<int>();
    p = new Point(dimension_ + 1);
    for(unsigned int j = 0; j < dimension_; ++j)
        (*p)[j] = 0;
    
    (*p)[dimension_ - 1] = initial_value[dimension_ - 1];
    (*p)[dimension_] = 1;
    
    for(unsigned int i = 0; i < dimension_ - 1; ++i)
        node_inequality_indices_[n]->push_back(i);
    
    node_inequality_indices_[n]->push_back(dimension_);
    birth_index_[n] = dimension_;
    
    point_nodes_.insert(make_pair(p, n));
    node_points_[n] = p;
    unprocessed_projective_points_.push(p);
    
    forall_nodes(v, vertex_graph_) {
        if(v != n) {
            vertex_graph_.newEdge(v, n);
            vertex_graph_.newEdge(n, v);
        }
    }
    
    p = new Point(dimension_ + 1);
    for(unsigned int j = 0; j < dimension_; ++j)
        (*p)[j] = initial_value[j] - initial_value[dimension_ - 1];
    (*p)[dimension_ - 1] = -1;
    (*p)[dimension_] = initial_value[dimension_ - 1];
    
    list_of_inequalities_.push_back(p);
    
    n = vertex_graph_.newNode();
    node_inequality_indices_[n] = new list<int>();
    p = new Point(dimension_ + 1);
    for(unsigned int j = 0; j < dimension_ - 1; ++j)
        (*p)[j] = 0;
    (*p)[dimension_ - 1] = -1;
    (*p)[dimension_] = 0;
    
    for(unsigned int i = 0; i < dimension_; ++i)
        node_inequality_indices_[n]->push_back(i);
    
    birth_index_[n] = dimension_ - 1;
    
    point_nodes_.insert(make_pair(p, n));
    node_points_[n] = p;
    
    forall_nodes(v, vertex_graph_) {
        if(v != n) {
            vertex_graph_.newEdge(v, n);
            vertex_graph_.newEdge(n, v);
        }
    }
    
    //	cout << "inequalities:" << endl;
    //	for(auto ineq : list_of_inequalities_)
    //		cout << *ineq << endl;
    //
    //	cout << "vertices:" << endl;
    //	forall_nodes(v, vertex_graph_) {
    //		cout << *node_points_[v] << endl;
    //		cout << "Node inequalities:" << endl;
    //			for(auto index: *node_inequality_indices_[v])
    //				cout << index << ", ";
    //			cout << endl;
    //	}
    //
    //	cout << "graph has " << vertex_graph_.numberOfNodes() << " nodes and " << vertex_graph_.numberOfEdges() << " edges" << endl;
}

NodeListVE::~NodeListVE() {
	node n;

	forall_nodes(n, vertex_graph_) {
		delete node_points_[n];
	}

	for(auto inequality : list_of_inequalities_)
		delete inequality;

}

inline Point* NodeListVE::
to_projective(Point &non_projective_point) {
    
	double * projective_values = new double[dimension_ + 1];

	for(unsigned int i = 0; i < dimension_; ++i)
		projective_values[i] = non_projective_point[i];

	projective_values[dimension_] = 1;

	Point *projective_point = new Point(projective_values, dimension_ + 1);
	delete[] projective_values;

	return projective_point;
}

inline Point NodeListVE::
normalize_projective(Point projective_point) {

	double *normalized_projective_values = new double[dimension_ + 1];

	for(unsigned int i = 0; i < dimension_; ++i)
		normalized_projective_values[i] = projective_point[i] / projective_point[dimension_];

	normalized_projective_values[dimension_] = 1;

	Point normalized_projective_point(normalized_projective_values, dimension_ + 1);
	delete[] normalized_projective_values;

	return normalized_projective_point;
}

inline node NodeListVE::
get_node(Point &non_projective_point) {
    
	Point *projective_point = to_projective(non_projective_point);
	node n = point_nodes_[projective_point];

	assert(EqualityPointComparator(epsilon_)(node_points_[n], projective_point));

	delete projective_point;


	return n;
}

Point * NodeListVE::next_vertex() {
	Point *point;

//	cout << *unprocessed_projective_points_.front() << " noch im Graphen: " << point_nodes_.count(unprocessed_projective_points_.front()) << endl;
	while(!unprocessed_projective_points_.empty() && point_nodes_.count(unprocessed_projective_points_.top()) == 0)
		unprocessed_projective_points_.pop();

//	cout << *node_points_[point_nodes_[unprocessed_projective_points_.front()]] << endl;

	if(unprocessed_projective_points_.empty())
		return nullptr;

	point = unprocessed_projective_points_.top();
//	cout << *point << ", dim: " << point->dimension() << endl;
	Point * new_point = new Point(dimension_);

//	for(auto entity: point_nodes_)
//		cout << *entity.first << " - " << entity.second << ", dim: " << entity.first->dimension() << endl;

	assert(point->dimension() == dimension_ + 1);

	if(abs((*point)[dimension_] - 1) > epsilon_) {
		cout << "Point " << *point << " is not normalized: " << (*point)[dimension_] - 1 << endl;
		assert(false);
	}

	for(unsigned int i = 0; i < dimension_; ++i)
		(*new_point)[i] = (*point)[i];

	unprocessed_projective_points_.pop();

	return new_point;
}

void NodeListVE::add_hyperplane(Point &vertex, Point &normal, double rhs) {
//	cout << vertex << endl;
//	cout << normal << endl;
//	cout << rhs << endl;

	clock_t start = clock();

	Point *projective_normal = to_projective(normal);
	(*projective_normal)[dimension_] = -rhs;

	list_of_inequalities_.push_back(projective_normal);
//	cout << *projective_normal << endl;

	// Check if redundancy is introduced
//	dd_set_global_constants();
//	dd_MatrixPtr A;
//	dd_ErrorType err;
//	dd_rowset rows;
//	A = dd_CreateMatrix(list_of_inequalities_.size(), dimension_ + 1);
//	int j = 0;
//	for(auto row: list_of_inequalities_) {
//
//		for(unsigned int i = 0; i < dimension_; ++i)
//			dd_set_d(A->matrix[j][i + 1], (*row)[i]);
//
//		dd_set_d(A->matrix[j][0], (*row)[dimension_]);
//
//		j+= 1;
//	}
//	A->representation = dd_Inequality;
//	rows = dd_RedundantRows(A, &err);
//
//	if(set_card(rows) > 0) {
//		cout << "CDD giving redundant inequalities: " << set_card(rows) << endl;
//		dd_WriteMatrix(stdout, A);
//		set_write(rows);
//		assert(false);
//	}
//	set_free(rows);
//	dd_FreeMatrix(A);
//	dd_free_global_constants();
	// End Redundancy check


	list<node> active_nodes;
	active_nodes.push_back(get_node(vertex));

	list<node> new_face_nodes;

	NodeArray<bool> already_active(vertex_graph_, false);
	NodeArray<node> father_node(vertex_graph_, nullptr);

	bool nondegenerate = true;

	while(!active_nodes.empty()) {

		node active_node = active_nodes.front();
		active_nodes.pop_front();
//		cout << "active node: " << active_node << endl;
		assert(active_node != nullptr);

		Point *projective_vertex = node_points_[active_node];

//		cout << "vertex: " << *projective_vertex << endl;

		NodeArray<bool> neighor_checked(vertex_graph_, false);

		AdjElement *adj;
		forall_adj(adj, active_node) {

			node neighbor = adj->theEdge()->target();

//			cout << "neighbor: " << neighbor << endl;

			if(neighbor == active_node)
				continue;

			if(neighor_checked[neighbor])
				continue;

			neighor_checked[neighbor] = true;

			Point *projective_check_point = node_points_[neighbor];

//			cout << "neighbor to check: " << *projective_check_point << endl;

//			for(auto n : new_face_nodes)
//				cout << "1 face node: " << n << " with point: " << *node_points_[n] << " at address " << (void *) n << endl;

			node new_node = nullptr;

			// Check if neighbor is on the other side of the hyperplane
			if((*projective_normal) * (*projective_check_point) < -epsilon_) {
				if(!already_active[neighbor]) {
					active_nodes.push_back(neighbor);
					already_active[neighbor] = true;
				}
//				cout << "neighbor is also outside" << endl;
				continue;

			} else if((*projective_normal) * (*projective_check_point) < epsilon_) {
				nondegenerate = false;
				new_node = neighbor;
//				cout << "neighbor " << new_node << " resides on hyperplane" << endl;
				bool already_on_face = false;
				for(auto n : new_face_nodes) {
					if(new_node->index() == n->index())
						already_on_face = true;
				}
				if(already_on_face)
					continue;

			} else {

//				cout << "On hyperplane check: " << (*projective_normal) * (*projective_check_point) << endl;

				double alpha = - (*projective_normal * *projective_vertex);
				Point diff_direction = (*projective_check_point) - (*projective_vertex);

				if(abs((*projective_normal) * diff_direction) < epsilon_) {
//					cout << "divisor: " << (*projective_normal) * diff_direction << endl;
					continue;
				} else
					alpha /=  (*projective_normal) * diff_direction;

//				for(auto n : new_face_nodes)
//					cout << "2 face node: " << n << " with point: " << *node_points_[n] << " at address " << (void *) n << endl;

//				cout << "diff direction: " << diff_direction << endl;
//				cout << "alpha: " << alpha << endl;

				assert( alpha > epsilon_ && alpha < 1 + epsilon_);

//				for(auto n : new_face_nodes)
//					cout << "3 face node: " << n << " with point: " << *node_points_[n] << " at address " << (void *) n << endl;

				diff_direction *= alpha;
				Point *projective_cut_point = new Point(normalize_projective((*projective_vertex) + diff_direction));

//				cout << "new cut vertex: " << *projective_cut_point << endl;

				unprocessed_projective_points_.push(projective_cut_point);

//				for(auto n : new_face_nodes)
//					cout << "4 face node: " << n << " with point: " << *node_points_[n] << " at address " << (void *) n << endl;

//				cout << adj->theEdge() << endl;
				edge new_edge = vertex_graph_.split(adj->theEdge());
//				cout << adj->theEdge() << endl;
				new_node = new_edge->source();
				node_inequality_indices_[new_node] = new list<int>();
//				cout << (void *) new_node << endl;

				set_intersection(	node_inequality_indices_[active_node]->begin(),
									node_inequality_indices_[active_node]->end(),
									node_inequality_indices_[neighbor]->begin(),
									node_inequality_indices_[neighbor]->end(),
									back_inserter(*node_inequality_indices_[new_node]));

//				cout << "New node: " << new_node << endl;

				vertex_graph_.newEdge(neighbor, new_node);

				node_points_[new_node] = projective_cut_point;
				father_node[new_node] = neighbor;
				birth_index_[new_node] = list_of_inequalities_.size() - 1;

//				for(auto n : new_face_nodes)
//					cout << "5 face node: " << n << " with point: " << *node_points_[n] << " at address " << (void *) n << endl;

				point_nodes_.insert(make_pair(projective_cut_point, new_node));
			}

//			for(auto n : new_face_nodes)
//				cout << "face node: " << n << " with point: " << *node_points_[n] << " at address " << (void *) n << endl;

			new_face_nodes.push_back(new_node);
			node_inequality_indices_[new_node]->push_back(list_of_inequalities_.size() - 1);

//			for(auto n : new_face_nodes)
//				cout << "face node: " << n << " with point: " << *node_points_[n] << " at address " << (void *) n << endl;

		}	// forall_adj

//		cout << "Deleting node: " << active_node << endl;
		delete node_inequality_indices_[active_node];
		vertex_graph_.delNode(active_node);
		point_nodes_.erase(projective_vertex);
//		delete projective_vertex;

	}	// active node

//	cout << "Checking face nodes: ";
//	for(auto n : new_face_nodes)
//		cout << n << ", ";
//	cout << endl;

	for(auto iter1 = new_face_nodes.begin(); iter1 != new_face_nodes.end(); ++iter1) {
		auto n1 = *iter1;

		for(auto iter2 = iter1; iter2 != new_face_nodes.end(); ++iter2) {

			auto n2 = *iter2;

			if(n1 == n2)
				continue;

			if(father_node[n1] != n2 && father_node[n2] != n1)
				if(inside_face(n1, n2, nondegenerate))
					continue;

			vertex_graph_.newEdge(n1, n2);
			vertex_graph_.newEdge(n2, n1);
		}

		if(dimension_ == 2)
			assert((n1->indeg() >= 2 && n1->outdeg() >= 2) || abs((*node_points_[n1])[dimension_]) < epsilon_);
		else
			if(((unsigned int) n1->indeg() < dimension_ || (unsigned int) n1->outdeg() < dimension_) && abs((*node_points_[n1])[dimension_]) > epsilon_) {
				cout << "node " << n1 << " with point " << *node_points_[n1] << " has indegree " << n1->indeg() << " and outdegree " << n1->outdeg() << endl;
				assert(false);
			}
	}

	cycles_ += clock() - start;
}

bool NodeListVE::inside_face(node n1, node n2, bool nondegenerate) {
//	cout << "p1: " << p1 << " (" << point_nodes_[&p1] << "), p2: " << p2 << " (" << point_nodes_[&p2] << ")" << endl;

	Point p1 = *node_points_[n1];
	Point p2 = *node_points_[n2];

	if(dimension_ == 2)
		return false;

	else if(dimension_ > 3 || !nondegenerate) {

		assert(n1 != n2);

		if(abs(p1[dimension_]) < epsilon_ && abs(p2[dimension_]) < epsilon_)
			return false;

		set<Point *> common_vertices;

		unsigned int i = 0;
		list<unsigned int> tight_inequalities;

		set_intersection(	node_inequality_indices_[n1]->begin(),
							node_inequality_indices_[n1]->end(),
							node_inequality_indices_[n2]->begin(),
							node_inequality_indices_[n2]->end(),
							back_inserter(tight_inequalities));

		unsigned int k = max(birth_index_[n1], birth_index_[n2]);
		bool nc2 = false;

		for(auto index: tight_inequalities)
			if(index >= k && index <= list_of_inequalities_.size() - 1) {
				nc2 = true;
				break;
			}

		// [FP96] NC2
		if(!nc2)
			return true;

		// [FP96] NC1
		if(tight_inequalities.size() < dimension_ - 2)
			return true;

		for(auto inequality_index : tight_inequalities) {

			auto inequality = list_of_inequalities_[inequality_index];

			set<Point *> new_vertices;

//			cout << "inequality: " << *inequality << endl;

			node n;
			forall_nodes(n, vertex_graph_) {
//					cout << "checking node " << n << " with point " << *node_points_[n] << ": " << (*inequality) * (*node_points_[n]) << endl;
				if(abs((*inequality) * (*node_points_[n])) < epsilon_)
					new_vertices.insert(node_points_[n]);
			}

			if(common_vertices.empty())
				common_vertices.insert(new_vertices.begin(), new_vertices.end());
			else {
				list<Point *> temp_points;
				set_intersection(common_vertices.begin(), common_vertices.end(), new_vertices.begin(), new_vertices.end(), back_inserter(temp_points));
//					cout << "intersection size: "<< temp_points.size() << endl;
				common_vertices.clear();
				common_vertices.insert(temp_points.begin(), temp_points.end());
			}

//				cout << "current number of common vertices: " << common_vertices.size() << endl;

			i += 1;
		}

//		cout << "Number of common vertices: " << common_vertices.size() << endl;

		if(common_vertices.size() == 2)
			return false;
		else
			return true;

	} else if(dimension_ == 3 || nondegenerate) {

		list<int> inequality_intersection;

		node n1 = point_nodes_[&p1];
		node n2 = point_nodes_[&p2];

		set_intersection(	node_inequality_indices_[n1]->begin(),
							node_inequality_indices_[n1]->end(),
							node_inequality_indices_[n2]->begin(),
							node_inequality_indices_[n2]->end(),
							back_inserter(inequality_intersection));

		unsigned int tight_inequalities = inequality_intersection.size();

		if(tight_inequalities == dimension_ - 1)
			return false;
		else if(tight_inequalities < dimension_ - 1)
			return true;
		else
			assert(false);
	}

	return false;
}

} /* namespace mco */
