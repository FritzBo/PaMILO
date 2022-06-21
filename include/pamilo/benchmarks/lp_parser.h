/**
 * @file  lp_parser.h
 * @author Mirko H. Wagner
 * @brief
 * @date 06.06.2020
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#pragma once

#include <string>

#include <pamilo/pilp/ilp.h>

namespace pamilo {

/**
 * @brief Class to read an ILP from a .lp file
 *
 */
class LPparser
{
public:
    /**
     * @brief Loads a .lp file, parses the ILP and stores it in ilp
     *
     * @param filename
     * @param ilp
     */
    void getILP(std::string filename, ILP &ilp);
};
}  // namespace pamilo
