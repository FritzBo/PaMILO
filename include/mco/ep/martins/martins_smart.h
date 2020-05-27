#pragma once
/*

 * martins.h
 *
 *  Created on: 12.03.2013
 *      Author: fritz
 */

#ifndef MARTINS_S_H_
#define MARTINS_S_H_

#include <mco/ep/basic/abstract_ep_solver.h>

namespace mco {

class EpSolverMartinsSmart : public AbstractEpSolver {

public:
	explicit EpSolverMartinsSmart(mco::EpInstance &instance);
	virtual void Solve();

};

}

#endif /* MARTINS_S_H_ */
