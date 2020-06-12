#pragma once
/*
 * graph_generators.h
 *
 *  Created on: 28.04.2013
 *      Author: fritz
 */

#ifndef GRAPH_GENERATORS_H_
#define GRAPH_GENERATORS_H_

#include <ogdf/basic/Graph.h>

namespace mco {

void KTreeGenerator(ogdf::Graph & graph, unsigned int number_of_nodes);

}

#endif /* GRAPH_GENERATORS_H_ */
