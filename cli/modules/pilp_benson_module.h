/**
 * @file pilp_benson_module.h
 * @author Fritz Bökler and Mirko H. Wagner
 * @brief
 * @date 28.05.2020
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#pragma once

#include <list>
#include <string>

#include "../basic/modules.h"

/**
 * @brief Module to calculate non-dominated extreme points of a MOMILP
 *
 */
class PilpBensonModule : public AlgorithmModule<std::string>
{
public:
    /**
     * @brief
     *
     * @param argc
     * @param args
     */
    virtual void perform(int argc, char **args);
    virtual ~PilpBensonModule()
    {
    }

    virtual const std::list<std::pair<const std::string, const pamilo::Point>> &solutions();
    virtual std::string statistics();

private:
    std::list<std::pair<const std::string, const pamilo::Point>> solutions_;
};
