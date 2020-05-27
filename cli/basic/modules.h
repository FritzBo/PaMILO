//
//  modules.h
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#ifndef __mco__modules__
#define __mco__modules__

#include <ostream>
#include <string>
#include <map>
#include <list>
#include <ostream>

#include <tclap/CmdLine.h>

#include <mco/basic/point.h>

class BasicModule {
public:
    virtual void perform(int argc, char** args) = 0;
    virtual ~BasicModule() = default;
};

template<class T>
class AlgorithmModule : public BasicModule {
public:
    
    using solution_type = T;
    using csolution_type = const T;
    using solution_type_pointer = T*;
    using csolution_type_pointer = const T*;
    
    virtual const std::list<std::pair<csolution_type, const mco::Point>>& solutions() = 0;
    virtual std::string statistics() = 0;
};

class ModuleFactory {
public:
    void add_module(std::string name, BasicModule& module);
    
    std::list<std::pair<unsigned, BasicModule*>> parse_module_list(int argc, char** argv);
    
    ~ModuleFactory();

private:
    std::map<std::string, BasicModule*> modules_;
};
#endif /* defined(__mco__modules__) */
