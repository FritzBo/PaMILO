//
//  pilp_benson_module.h
//  pamilo
//
//  Created by Fritz Bökler and Mirko H. Wagner on 28.05.20.
//
//

#ifndef __pamilo__pilp_benson_module__
#define __pamilo__pilp_benson_module__

#include <list>

#include "../basic/modules.h"

class PilpBensonModule : public AlgorithmModule<std::string> {

public:
    virtual void perform(int argc, char** args);
    virtual ~PilpBensonModule() {}

    virtual const std::list<std::pair<const std::string, const pamilo::Point>>& solutions();
    virtual std::string statistics();

private:

    std::list<std::pair<const std::string, const pamilo::Point>> solutions_;


};

#endif /* defined(__pamilo__pilp_benson_module__) */
