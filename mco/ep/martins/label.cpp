/*
 * label.cpp
 *
 *  Created on: 27.03.2013
 *      Author: fritz
 */

#include<mco/ep/martins/label.h>

namespace mco {

Label::Label(const Point *point, ogdf::node n, const Label *pred) : point(point), n(n), pred(pred), mark_dominated(false) {

}

Label::Label(const Label &label) : point(new Point(*label.point)), n(label.n), pred(label.pred), mark_dominated(label.mark_dominated) {

}

Label::~Label() {
	delete point;
}

}

