//
//  pilp_benson_module.h
//  mco
//
//  Created by Fritz Bökler and Mirko H. Wagner on 28.05.20.
//
//

#ifndef __mco__pilp_benson_module__
#define __mco__pilp_benson_module__

#include <list>

#include "../basic/modules.h"

class PilpBensonModule : public AlgorithmModule<std::list<std::string>> {
    
public:
    virtual void perform(int argc, char** args);
    virtual ~PilpBensonModule() {}
    
    virtual const std::list<std::pair<const std::list<std::string>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:
    
    std::list<std::pair<const std::list<std::string>, const mco::Point>> solutions_;

    
};

#endif /* defined(__mco__pilp_benson_module__) */
