/*
 * correlated_random.cpp
 *
 *  Created on: 12.04.2013
 *      Author: fritz
 */

#include <ogdf/basic/Graph.h>

using ogdf::randomDouble;

#include <mco/benchmarks/correlated_random.h>
#include <mco/basic/point.h>

namespace mco {

Point * CorrelatedObjectiveGenerator::draw_point() {
	double d = randomDouble(lower_bound_, upper_bound_);
	double d2;
	if(rho_ >= 0)
		d2 = (rho_ * d) + (1 - rho_) * randomDouble(lower_bound_, upper_bound_);
	else {
		d2 = (upper_bound_ - d) * -rho_ + (1 + rho_) * randomDouble(lower_bound_, upper_bound_);
	}
	double *values = new double[2];
	values[0] = d;
	values[1] = d2;
	return new Point(values, 2);
}

}
