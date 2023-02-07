/**
 * @file interface_util.hpp
 * @author Levin Nemesch
 * @brief
 * @date 28.05.2020
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */
#pragma once

namespace pamilo
{
    enum SolverStatus{Unbounded, Infeasible, Success, Error, UnknownStatus};

    enum VarType{Float, Integer, UnknownType};
} // namespace pamilo
