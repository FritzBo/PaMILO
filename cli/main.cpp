/**
 * @file main.cpp
 * @author Fritz Bökler
 * @brief
 * @date 09.04.2014
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#include <iostream>

using std::cout;
using std::endl;
using std::exception;

#include <pamilo/basic/point.h>
#include <pamilo/pilp/pilp_dual_benson.h>

using pamilo::ComponentwisePointComparator;
using pamilo::EqualityPointComparator;
using pamilo::LexPointComparator;
using pamilo::ParetoDominationPointComparator;
using pamilo::Point;

#include "basic/modules.h"
#include "modules/pilp_benson_module.h"

int main(int argc, char **argv)
{
    try
    {
        PilpBensonModule pilp_benson_module;

        pilp_benson_module.perform(argc, argv);
    }
    catch (exception &e)
    {
        cout << e.what() << endl;
    }
}
