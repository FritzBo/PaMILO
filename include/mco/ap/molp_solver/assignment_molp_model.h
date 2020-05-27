#pragma once
/*
 * assignment_molp_solver.h
 *
 *  Created on: 19.08.2013
 *      Author: fritz
 */

#ifndef ASSIGNMENT_MOLP_SOLVER_H_
#define ASSIGNMENT_MOLP_SOLVER_H_

#include <memory>

#include <mco/molp/basic/molp_model.h>
#include <mco/ap/basic/ap_instance.h>

namespace mco {

class AssignmentMolpModel : public MolpModel {
public:
	AssignmentMolpModel(std::shared_ptr<AssignmentInstance> instance);

private:
	std::shared_ptr<AssignmentInstance> instance_;
};

} /* namespace mco */
#endif /* ASSIGNMENT_MOLP_SOLVER_H_ */
