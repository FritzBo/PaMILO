/*
 * random_objective_generator.cpp
 *
 *  Created on: 12.08.2013
 *      Author: fritz
 */

#include <ogdf/basic/Graph.h>

using ogdf::randomDouble;

#include <mco/basic/point.h>
#include <mco/benchmarks/random_objective_generator.h>

namespace mco {

Point * RandomObjectiveGenerator::draw_point() {
	double *values = new double[dimension()];
	for(unsigned int i = 0; i < dimension(); ++i)
		values[i] = randomDouble(lower_bound_, upper_bound_);

	Point *point = new Point(values, dimension());
	delete[] values;

	points_.push_back(point);
	return point;
}

RandomObjectiveGenerator::~RandomObjectiveGenerator() {

}

} /* namespace mco */
