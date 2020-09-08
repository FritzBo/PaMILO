//
//  pilp_benson_module.h
//  pamilo
//
//  Created by Fritz Bökler and Mirko H. Wagner on 28.05.20.
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
