/*
 * graph_vertex_container.cpp
 *
 *  Created on: 14.08.2013
 *      Author: fritz
 */

#include <mco/molp/benson/graph_vertex_container.h>

#include <utility>
#include <cassert>
#include <cmath>
#include <set>
#include <algorithm>

using std::make_pair;
using std::list;
using std::abs;
using std::set;
using std::set_intersection;
using std::back_inserter;

using ogdf::node;
using ogdf::edge;
using ogdf::AdjElement;

using ogdf::NodeArray;

#include <setoper.h>
#include <cdd.h>

#include <mco/core/point.h>

namespace mco {

GraphVertexContainer::GraphVertexContainer(Point & ideal_point,
                                           unsigned int dimension,
                                           double epsilon)
    :   OnlineVertexEnumerator(dimension, epsilon) {
        
	node_points_.init(vertex_graph_);

	node ideal_node = vertex_graph_.newNode();
	Point * projective_ideal_point = to_projective(ideal_point);

	node_points_[ideal_node] = projective_ideal_point;
	point_nodes_.insert(make_pair(projective_ideal_point, ideal_node));

	Point *infinity_point, *initial_inequality;
	for(unsigned int i = 0; i < dimension; ++i) {
		double *values = new double[dimension + 1];

		for(unsigned int j = 0; j < dimension + 1; ++j)
			if(i == j)
				values[j] = 1;
			else
				values[j] = 0;

		infinity_point = new Point(values, dimension + 1);

		node n = vertex_graph_.newNode();
		node_points_[n] = infinity_point;
		point_nodes_.insert(make_pair(infinity_point, n));

		vertex_graph_.newEdge(ideal_node, n);

		values[dimension] = -ideal_point[i];
		initial_inequality = new Point(values, dimension + 1);
		list_of_inequalities_.push_back(initial_inequality);

		delete[] values;
	}

	unprocessed_projective_points_.push(projective_ideal_point);

}

} /* namespace mco */
