#pragma once
/*
 * pge10_parser3_ap.h
 *
 *  Created on: 06.11.2013
 *      Author: fritz
 */

#ifndef PGE10_PARSER3_AP_H_
#define PGE10_PARSER3_AP_H_

#include <memory>

#include <mco/ap/basic/ap_instance.h>

namespace mco {

class PGE10Parser3AP {
public:
	PGE10Parser3AP() = delete;

	PGE10Parser3AP(std::string filename) :
		filename_(filename) {}

	AssignmentInstance * get_instance();

private:
	std::string filename_;
};

} /* namespace mco */
#endif /* PGE10_PARSER3_AP_H_ */
