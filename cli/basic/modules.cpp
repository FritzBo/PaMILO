/**
 * @file modules.cpp
 * @author Fritz Bökler
 * @brief
 * @date 14.07.2014
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#include "modules.h"

using std::list;
using std::make_pair;
using std::pair;
using std::string;

void ModuleFactory::add_module(string name, BasicModule &module)
{
    modules_.insert(make_pair(name, &module));
}

ModuleFactory::~ModuleFactory()
{
}

list<pair<unsigned, BasicModule *>> ModuleFactory::parse_module_list(int argc, char **argv)
{
    list<pair<unsigned, BasicModule *>> found_modules;

    for (unsigned i = 0; i < argc; ++i)
    {
        auto module_it = modules_.find(string(argv[i]));

        if (module_it != modules_.end())
        {
            found_modules.push_back(make_pair(i, module_it->second));
        }
    }

    return found_modules;
}
