#pragma once
/*
 * abstract_vertex_container.h
 *
 *  Created on: 14.08.2013
 *      Author: fritz
 */

#ifndef ABSTRACT_VERTEX_CONTAINER_H_
#define ABSTRACT_VERTEX_CONTAINER_H_

#include <mco/point.h>

#include <list>

namespace mco {

class AbstractVertexContainer {
public:
	AbstractVertexContainer() {}

	AbstractVertexContainer(AbstractVertexContainer *) = delete;
	AbstractVertexContainer * operator=(AbstractVertexContainer *) = delete;

	virtual ~AbstractVertexContainer() {};
	virtual bool has_next() = 0;
	virtual Point * next_vertex() = 0;
	virtual void add_hyperplane(Point &vertex, Point &normal, double lhs) = 0;

	virtual unsigned int number_of_hyperplanes() = 0;
};

} /* namespace mco */
#endif /* ABSTRACT_VERTEX_CONTAINER_H_ */
