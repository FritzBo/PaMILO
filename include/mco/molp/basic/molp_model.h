#pragma once
/*
 * molp_model.h
 *
 *  Created on: 20.06.2013
 *      Author: fritz
 */

#ifndef MOLP_MODEL_H_
#define MOLP_MODEL_H_

#include <valarray>
#include <vector>

namespace mco {

class MolpModel {

public:
	MolpModel() :
		A_(nullptr), b_(nullptr), C_(nullptr), objectives_(0), constraints_(0), variables_(0) { }

	MolpModel(double *A, double *b, double *C, unsigned int objectives, unsigned int constraints, unsigned int variables) :
		A_(A), b_(b), C_(C), objectives_(objectives), constraints_(constraints), variables_(variables) { }

	const double *A() const {
		return A_;
	}

	const double *b() const {
		return b_;
	}

	const double *C() const {
		return C_;
	}

	unsigned int objectives() const {
		return objectives_;
	}

	unsigned int constraints() const {
		return constraints_;
	}

	unsigned int variables() const {
		return variables_;
	}

	friend std::ostream & operator<<(std::ostream &, const MolpModel &);

protected:
	double *A_;
	double *b_;
	double *C_;

	unsigned int objectives_;
	unsigned int constraints_;
	unsigned int variables_;


};

std::ostream & operator<<(std::ostream &os, const MolpModel &model);

}

#endif /* MOLP_MODEL_H_ */
