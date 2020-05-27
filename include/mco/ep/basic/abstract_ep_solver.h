#pragma once
/*
 * EpSolver.h
 *
 *  Created on: 15.03.2013
 *      Author: fritz
 */

#ifndef EPSOLVER_H_
#define EPSOLVER_H_

#include <list>

#include <ogdf/basic/Graph.h>

#include <mco/basic/abstract_solver.h>

namespace mco {

class Point;
class EpInstance;

class AbstractEpSolver : public AbstractSolver<std::list<ogdf::edge>> {

	EpInstance &instance_;

protected:

	EpInstance & instance() const { return instance_; }

public:
	AbstractEpSolver() = delete;
	explicit AbstractEpSolver(EpInstance &instance) : instance_(instance) {}

};

}

#endif /* EPSOLVER_H_ */
