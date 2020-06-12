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
		ifstream orig(filename);
		char temp[] = ".temp.lp";
		ofstream woObj(temp);
		string line;
		vector<string> objs;

		if(!orig.is_open()) {
			cerr << "input file not openable!\n";
		}
		if(!woObj.is_open()) {
			cerr << "cannot create tempfile\n";
		}

		while(getline(orig, line)) {
			cout << line << endl;
			// TODO: maybe case insensitive 
			if(strstr(line.c_str(), "aximize") || strstr(line.c_str(), "inimize")) {
				// were in the obj function section now!
				woObj << line << endl;
				// a dummy objective is pushed, so it's not empty
				woObj << "obj:x1\n";
				int dimension = 0;
				while(getline(orig, line) && strstr(line.c_str(), "Subject To")) {
					// a line containing a name followed by an ':' indicates a new objective
					if(strstr(line.c_str(), ":")){
						++dimension;
					}
					objs.push_back(line);
				}
			} else {
				woObj << line << endl;
			}
		}
		orig.close();
		woObj.close();


		// read new file
		ilp.osi.readLp(temp);

		// create map to get column index by name
		map<string, int> mapColNameToIndex;
		for(int i = 0; i < ilp.osi.getNumCols(); i++) {
			mapColNameToIndex.insert(pair<string, int>(ilp.osi.getColName(i), i));
		}

		// read objectives
		int curDim = 0;
		for(string objLine : objs) {
			// read name of objective followed by a ':'
			if(!strstr(objLine.c_str(), ":")) {
				// TODO: split objLine into pairs of name and coeff
				// while {
					string name = "dummy";
					double coef = 0.0;
					int index = mapColNameToIndex[name];
					ilp.obj[curDim].insert(pair<int, double>(index, coef));
				// }
				curDim++;
			}
			// we assume no parameters were given in a new line (in the name line it's ok)
		}


	}
}
