//
//  lp_parser.h
//  pamilo
//
//  Created by Mirko H. Wagner on 06.06.20.
//
//  This file is distributed for academics only
//  under the terms of an MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
//
//

#pragma once

#include <string>

#include <pamilo/pilp/ilp.h>

namespace pamilo {

class LPparser {
public:
    void getILP(std::string filename,
                  ILP &ilp);
};
}

