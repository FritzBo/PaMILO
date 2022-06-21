/**
 * @file abstract_solver.h
 * @author Fritz Bökler
 * @brief
 * @date 25.03.2013
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#pragma once

#include <iterator>
#include <list>

#include <pamilo/basic/point.h>

namespace pamilo {

/**
 * @brief Solver base class. Internally stores a set of (found) solutions, provides no solving
 * capabilities.
 *
 * @tparam T Type in which solutions are represented
 */
template <class T>
class AbstractSolver
{
public:
    using solution_type = T;
    using csolution_type = const T;

    AbstractSolver() = default;
    AbstractSolver(const AbstractSolver &) = delete;

    AbstractSolver &operator=(const AbstractSolver &) = delete;

    /**
     * @brief Allows access to the stored set of solutions
     *
     * @return const std::list<std::pair<csolution_type, const Point>>&
     */
    const std::list<std::pair<csolution_type, const Point>> &solutions() const
    {
        return solutions_;
    }

    virtual ~AbstractSolver() = default;

protected:
    /**
     * @brief Adds a solution to the stored set of solutions
     *
     * @param solution
     * @param value
     */
    void add_solution(csolution_type solution, const Point value)
    {
        solutions_.push_back(std::make_pair(solution, value));
    }

    /**
     * @brief Add a number of solutions to the stored set of solutions
     *
     * @tparam InputIterator
     * @param begin
     * @param end
     */
    template <class InputIterator>
    void add_solutions(InputIterator begin, InputIterator end)
    {
        solutions_.insert(solutions_.end(), begin, end);
    }

    /**
     * @brief Clears the stored set of solutions
     *
     */
    void reset_solutions()
    {
        solutions_.clear();
    }

private:
    std::list<std::pair<csolution_type, const Point>> solutions_;
};
}  // namespace pamilo
