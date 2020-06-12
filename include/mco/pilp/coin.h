//
//  pilp.h 
//  mco
// 
//  Add
//
//  Created by Mirko H. Wagner on 20.06.20
//
//

#pragma once

#include <ogdf/basic/List.h>

// Cplex
#include <coin/OsiCpxSolverInterface.hpp>
typedef OsiCpxSolverInterface OSI;

// Glpk
//#include <coin/OsiGlpkSolverInterface.hpp>
//typedef OsiGlpkSolverInterface OSI;


class ILP {
public:
	OSI osi;
	int dimension;
	std::vector<std::map<int,double>> obj;

	ILP() {}
};
