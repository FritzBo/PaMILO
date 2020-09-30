//
//  lp_parser.h
//  pamilo
//
//  Created by Mirko H. Wagner on 06.06.20.
//
//  This file is distributed under the terms of
//
//  the GNU General Public License v3,
//  a copy of which can be found in the file LICENCE-GPLv3.txt
//
//  OR
//
//  for academics, a MIT license based license,
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

#endif /* defined(__pamilo__lp_parser__) */
