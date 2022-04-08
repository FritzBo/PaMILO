/**
 * @file lp_parser.cpp
 * @author Mirko H. Wagner
 * @brief
 * @date 20.06.2020
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#include <pamilo/benchmarks/lp_parser.h>

#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

namespace pamilo {
void LPparser::getILP(string filename, ILP &ilp)
{
    try
    {
        // ilp.cplex.importModel(ilp.model, filename.c_str(), ilp.obj, ilp.vars, ilp.cons);
        ilp.model = std::make_unique<GRBModel>(&ilp.env, filename);
        ilp.dimension = ilp.model->get(GRB_IntAttr_NumObj);
        ilp.obj = std::vector<GRBLinExpr>(ilp.dimension);
        ilp.multiObj = std::vector<GRBLinExpr>(ilp.dimension);
        for(int i=0; i<ilp.dimension; i++)
        {
            ilp.obj[i] = ilp.model->getObjective(i);
            ilp.multiObj[i] = ilp.model->getObjective(i);
        }
        ilp.vars = ilp.model->getVars();
        ilp.n_vars = ilp.model->get(GRB_IntAttr_NumVars);
    }
    catch (GRBException &e)
    {
        cerr << "GUROBI failed to read the file with error " << e.getMessage() << std::endl;;
        exit(-1);
    }
    // ilp.cplex.setParam(IloCplex::Param::MultiObjective::Display, 2);
    // ilp.cplex.setParam(IloCplex::Param::ParamDisplay, 0);


    ilp.relScale.resize(ilp.dimension, 1);
    ilp.offset.resize(ilp.dimension, 0);

    //ilp.cplex.extract(ilp.model);

    ilp.filename = filename;
}
}  // namespace pamilo
