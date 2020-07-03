//
//  lp_parser.h
//  mco
//
//  Created by Mirko H. Wagner on 06.06.20.
//
//

#ifndef __mco__lp_parser__
#define __mco__lp_parser__

#include <string>

#include <mco/pilp/coin.h>

namespace mco {
    
class LPparser {
public:
    void getILP(std::string filename,
                  ILP &ilp);
};
    
}

#endif /* defined(__mco__lp_parser__) */
