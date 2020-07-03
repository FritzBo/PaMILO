//
//  lp_parser.cpp
//  mco
//
//  Created by Mirko H. Wagner on 20.06.20
//
//

#include <mco/benchmarks/lp_parser.h>

#include <string.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

namespace mco {
	void LPparser::getILP(string filename, ILP &ilp) {
		// strip objectives from file (copy w/o obj into new file)
		ilp.osi.setIntParam(OsiIntParam::OsiNameDiscipline, 2);
		ilp.osi.readLp(filename.c_str());

		int objNum = 0;
		const double *objCols = ilp.osi.getObjCoefficients();
		for(int i = 0; i < ilp.osi.getNumCols(); i++) {
			if(objCols[i] != 0) {
				ilp.obj.push_back(i);
				objNum++;
				cout << ilp.osi.getColName(i) << endl;
			}
		}
		ilp.dimension = objNum;
		//cout << objNum << endl;
	}
}
