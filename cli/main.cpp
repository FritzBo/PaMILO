//
//  main.cpp
//  cli
//
//  Created by Fritz Bökler on 09.04.14.
//
//

#include <iostream>

using std::cout;
using std::endl;
using std::exception;

#include <mco/basic/point.h>
#include <mco/pilp/pilp_dual_benson.h>
#include <mco/generic/benson_dual/ove_cdd.h>

using mco::Point;
using mco::LexPointComparator;
using mco::ParetoDominationPointComparator;
using mco::ComponentwisePointComparator;
using mco::EqualityPointComparator;

#include "basic/modules.h"
#include "modules/pilp_benson_module.h"

int main(int argc, char** argv) {
	try {
		PilpBensonModule pilp_benson_module;
		
		pilp_benson_module.perform(argc, argv);
		
	} catch (exception& e) {
		cout << e.what() << endl;
	}
}

