//
//  modules.cpp
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#include "modules.h"

using std::string;
using std::list;
using std::make_pair;
using std::pair;

void ModuleFactory::add_module(string name, BasicModule& module) {
    modules_.insert(make_pair(name, &module));
}

ModuleFactory::~ModuleFactory() {
}

list<pair<unsigned, BasicModule*>> ModuleFactory::parse_module_list(int argc, char** argv) {
    list<pair<unsigned, BasicModule*>> found_modules;
    
    for(unsigned i = 0; i < argc; ++i) {
        auto module_it = modules_.find(string(argv[i]));
        
        if(module_it != modules_.end()) {
            found_modules.push_back(make_pair(i, module_it->second));
        }
    }
    
    return found_modules;
}