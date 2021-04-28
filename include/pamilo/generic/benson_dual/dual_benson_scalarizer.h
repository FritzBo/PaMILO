//
//  dual_benson_scalarizer.h
//
//  Created on: 27.09.2013
//      Author: Fritz Bökler
//
//  This file is distributed for academics only
//  under the terms of an MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
//
//

//#define WEPS
//#define VALEPS
//#define DMEAS

#pragma once

#include <list>
#include <functional>

#include <pamilo/basic/point.h>

namespace pamilo {

/*  OnlineVertexEnumerator Interface
 *
 *  OnlineVertexEnumerator(const Point& initial_value, unsigned dimension, double epsilon);
 *	bool has_next();
 *	Point * next_vertex();
 *	void add_hyperplane(Point &vertex, Point &normal, double rhs);
 *	unsigned int number_of_hyperplanes();
 */


template<typename OnlineVertexEnumerator, typename SolType = std::string>
class DualBensonScalarizer {
public:
	DualBensonScalarizer(
			std::function<double(const Point& weighting, Point& value, SolType &sol)> solver,
			std::function<void(std::pair<SolType, Point*>)> printSol,
			unsigned int dimension,
			double epsilon,
			double veEpsilon)
    :   dimension_(dimension),
		epsilon_(epsilon),
		veEpsilon_(veEpsilon),
		solver_(solver),
		printSol_(printSol),
		vertex_container(nullptr),
		vertices_(0),
		facets_(0),
		ve_time(0),
		oldWouldntButNewWould(0)
	{ }

	~DualBensonScalarizer() {
		if(vertex_container) {
			delete vertex_container;
		}
	}

	void Calculate_solutions(std::list<std::pair<SolType,Point *>>& solutions);

	double vertex_enumeration_time();

	int number_vertices() { return vertices_; }
	int number_facets() { return facets_; }

	int oldWouldntButNewWould;

protected:
	unsigned int dimension_;
	double epsilon_;
	double veEpsilon_;

	std::function<double(const Point&, Point&, SolType &sol)> solver_;
	std::function<void(std::pair<SolType, Point*>)> printSol_;

private:
	OnlineVertexEnumerator *vertex_container;

	int vertices_;
	int facets_;

	double ve_time;
};

template<typename OnlineVertexEnumerator, typename SolType>
void DualBensonScalarizer<OnlineVertexEnumerator, SolType>::
Calculate_solutions(std::list<std::pair<SolType, Point *>>& solutions) {
    vertices_ = 1;
    int iteration_counter = 0;
    facets_ = 1;

    Point v(dimension_);
    Point value(dimension_);

    for(unsigned int i = 0; i < dimension_ - 1; ++i)
        v[i] = 0;
	v[0] = 1;

	SolType sol;
	solver_(v, value, sol);
	auto solPair = std::make_pair(sol, new Point(value));
	solutions.push_back(solPair);
	printSol_(solPair);

	clock_t start = clock();
	if(vertex_container) {
		delete vertex_container;
	}
	vertex_container = new OnlineVertexEnumerator(value, dimension_, veEpsilon_);
	delete vertex_container->next_vertex();
	ve_time += (clock() - start) / (double) CLOCKS_PER_SEC;


	Point *candidate;
	Point weighting(dimension_);
	Point inequality(dimension_);

	double scalar_value;
    while(vertex_container->has_next()) {
        iteration_counter++;

#ifndef NDEBUG
        std::cout << "Iteration: " << iteration_counter << std::endl;
#endif

        candidate = vertex_container->next_vertex();

#ifndef NDEBUG
        std::cout << "New candidate: " << *candidate << std::endl;
#endif

        double sum = 0;
        for(unsigned int i = 0; i < dimension_ - 1; ++i) {
            weighting[i] = (*candidate)[i];
            sum += (*candidate)[i];
        }
        weighting[dimension_ - 1] = 1 - sum;

        for(unsigned int i = 0; i < dimension_; ++i) {
			value[i] = 0;
			if(weighting[i] < 0) {
#ifndef NDEBUG
				std::printf("\n\n\nweight eps: %E\n\n\n", weighting[i]);
#endif
#ifdef WEPS
				weighting[i] = 0;
#endif
			}
		}

#ifndef NDEBUG
        std::cout << "weighting: " << weighting << std::endl;
#endif

		SolType sol;
        scalar_value = solver_(weighting, value, sol);

#ifndef NDEBUG
        std::cout << "scalar value: " << scalar_value << std::endl;
        std::cout << "value vector: " << value << std::endl;
		std::cout << "candidate value: " << (*candidate)[dimension_ -1] << std::endl;
#endif
		for(unsigned int i = 0; i < dimension_ - 1; ++i) {
			double val = value[i] - value[dimension_ -1];
			if(abs(val) < epsilon_) {
#ifndef NDEBUG
				std::cout << "\n\n\nval eps. val: " << val << ", eps: " << epsilon_ << " component: " << i << ", from value[i] = " << value[i] << " and value[dim-1] = " << value[dimension_ -1] << "\n\n\n";
#endif
#ifdef VALEPS
				val = 0;
#endif
			}
			inequality[i] = val;
		}
		inequality[dimension_ - 1] = -1;

		// theoretically both values should be the same
		// first check is for debug, second is the real deal
		if(scalar_value - (*candidate)[dimension_ - 1] <= -epsilon_
				&& vertex_container->getDistance(*candidate,
				                                 inequality,
				                                 -value[dimension_ -1])> -epsilon_)
		{
			oldWouldntButNewWould++;
			std::cout << scalar_value - (*candidate)[dimension_ - 1] << " vs. "
				<< vertex_container->getDistance(*candidate, inequality, -value[dimension_ -1]) << std::endl;
		}
        if(scalar_value - (*candidate)[dimension_ - 1] > -epsilon_
#ifdef DMEAS
				|| vertex_container->getDistance(*candidate,
				                                 inequality,
												 -value[dimension_ -1]) > -epsilon_
#endif
				) {
			facets_++;
#ifndef NDEBUG
            std::cout << "found a new permanent extreme point. continuing." << std::endl;
#endif
        } else {
			clock_t start = clock();
            vertex_container->add_hyperplane(*candidate, inequality, -value[dimension_ - 1]);
			ve_time += (clock() - start) / (double) CLOCKS_PER_SEC;
			vertices_++;

			auto solPair = std::make_pair(sol, new Point(value));
			solutions.push_back(solPair);
			printSol_(solPair);
        }

        delete candidate;
    }

#ifndef NDEBUG
		std::cout << "Found " << facets_ << " nondominated value vectors in " << iteration_counter << " iterations." << std::endl;
		std::cout << "Where " << vertices_ << " weightings have been explored." << std::endl;
#endif
}

template<typename OnlineVertexEnumerator, typename SolType>
double DualBensonScalarizer<OnlineVertexEnumerator, SolType>::
vertex_enumeration_time() {
    return ve_time;
}
}

