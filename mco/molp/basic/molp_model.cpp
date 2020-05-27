/*
 * molp_model.cpp
 *
 *  Created on: 20.06.2013
 *      Author: fritz
 */

#include <mco/molp/basic/molp_model.h>

#include <ostream>

using std::ostream;
using std::endl;

namespace mco {

ostream & operator<<(ostream &os, const MolpModel &model) {
	os << "C:";

	for(unsigned int i = 0; i < model.objectives_; ++i) {
		os << "\t";
		for(unsigned int j = 0; j < model.variables_; ++j) {
			os << model.C_[i * model.variables_ + j] << ", ";
		}
		os << endl;
	}

	os << endl << "A:";
	for(unsigned int i = 0; i < model.constraints_; ++i) {
			os << "\t";
			for(unsigned int j = 0; j < model.variables_; ++j) {
				os << model.A_[i * model.variables_ + j] << ",\t";
			}
			os << endl;
		}

	os << endl << "b:\t";
	for(unsigned int i = 0; i < model.constraints_; ++i)
		os << model.b_[i] << ", ";

	return os;
}


}
