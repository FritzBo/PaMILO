#pragma once
/*
 * label.h
 *
 *  Created on: 27.03.2013
 *      Author: fritz
 */

#ifndef LABEL_H_
#define LABEL_H_

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

struct Label {
	const Point * const point;
	ogdf::node n;
	const Label * const pred;
	bool mark_dominated;

	Label(const Point *point, ogdf::node n, const Label *pred);
	Label(const Label &label);
	Label & operator=(const Label &label) = delete;

	~Label();
};

}

#endif /* LABEL_H_ */
