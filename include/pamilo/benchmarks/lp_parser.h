//
//  lp_parser.h
//  pamilo
//
//  Created by Mirko H. Wagner on 06.06.20.
//
//

#ifndef __pamilo__lp_parser__
#define __pamilo__lp_parser__

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
