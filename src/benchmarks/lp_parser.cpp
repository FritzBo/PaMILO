/**
 * @file lp_parser.cpp
 * @author Mirko H. Wagner and Levin Nemesch
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

#ifdef USE_GRB

void LPparser::getILP(string filename, ILP &ilp)
{
    try
    {
        ilp.model = std::make_unique<GRBModel>(&ilp.env, filename);
        ilp.dimension = ilp.model->get(GRB_IntAttr_NumObj);

        ilp.vars = std::unique_ptr<GRBVar[]>(ilp.model->getVars());
        ilp.n_vars = ilp.model->get(GRB_IntAttr_NumVars);

        // Overwrite tolerances from lp file to 0 (which is gurobi's default value here)
        for (int i = 0; i < ilp.dimension; i++)
        {
            ilp.model->set(GRB_DoubleAttr_ObjNRelTol, 0);
            ilp.model->set(GRB_DoubleAttr_ObjNAbsTol, 0);
        }

        // Always minimizing internally
        ilp.sense_og = ilp.model->get(GRB_IntAttr_ModelSense);
        if (ilp.sense_og == GRB_MAXIMIZE)
        {
            for (int i = 0; i < ilp.dimension; i++)
            {
                ilp.model->setObjectiveN(ilp.sense_og * ilp.model->getObjective(i), i);
            }
            ilp.model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        }

        ilp.model->set(GRB_StringParam_LogFile, ilp.grbFileName);
        ilp.model->update();
    }
    catch (GRBException &e)
    {
        cerr << "GUROBI failed to read the file with error " << e.getMessage() << std::endl;
        ;
        exit(-1);
    }

    ilp.relScale.resize(ilp.dimension, 1);
    ilp.offset.resize(ilp.dimension, 0);

    ilp.filename = filename;
}

#elif USE_CPLEX

void LPparser::getILP(string filename, ILP &ilp)
{
    try
    {
        ilp.cplex.importModel(ilp.model, filename.c_str(), ilp.obj, ilp.vars, ilp.cons);
        ilp.multiObj = ilp.obj;
    }
    catch (IloException &e)
    {
        cerr << "CPLEX failed to read the file. This is likely, because the input file is "
                "corrupted or is not a MOMIP\n";
        exit(-1);
    }
    ilp.cplex.setParam(IloCplex::Param::MultiObjective::Display, 2);
    ilp.cplex.setParam(IloCplex::Param::ParamDisplay, 0);
    ilp.cplex.setParam(IloCplex::Param::Threads, 1); // 1 as default
    ilp.cplex.setOut(ilp.cplexFile);

    ilp.dimension = ilp.multiObj.getNumCriteria();

    ilp.relScale.resize(ilp.dimension, 1);
    ilp.offset.resize(ilp.dimension, 0);

    ilp.cplex.extract(ilp.model);
    
    ilp.filename = filename;
}

#endif

}  // namespace pamilo
