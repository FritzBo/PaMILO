/**
 * @file modules.h
 * @author Fritz Bökler
 * @brief
 * @date 14.07.2014
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#pragma once

#include <list>
#include <map>
#include <ostream>
#include <string>

#include <tclap/CmdLine.h>

#include <pamilo/basic/point.h>

/**
 * @brief Base class for modules
 *
 */
class BasicModule
{
public:
    /**
     * @brief Performs the algorithm of the module
     *
     * @param argc
     * @param args
     */
    virtual void perform(int argc, char **args) = 0;
    virtual ~BasicModule() = default;
};

/**
 * @brief Base class for algorithmic solver modules
 *
 * @tparam T
 */
template <class T>
class AlgorithmModule : public BasicModule
{
public:
    using solution_type = T;
    using csolution_type = const T;
    using solution_type_pointer = T *;
    using csolution_type_pointer = const T *;

    /**
     * @brief Returns a list of solutions after they are calculated by perform()
     *
     * @return const std::list<std::pair<csolution_type, const pamilo::Point>>&
     */
    virtual const std::list<std::pair<csolution_type, const pamilo::Point>> &solutions() = 0;
    virtual std::string statistics() = 0;
};

class ModuleFactory
{
public:
    void add_module(std::string name, BasicModule &module);

    std::list<std::pair<unsigned, BasicModule *>> parse_module_list(int argc, char **argv);

    ~ModuleFactory();

private:
    std::map<std::string, BasicModule *> modules_;
};
