/**
 * @file ove_cdd.h
 * @author Fritz Bökler
 * @brief
 * @date 30.09.2013
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#pragma once

#include <iostream>
#include <list>

#include <cdd.h>
#include <setoper.h>

#include <pamilo/basic/point.h>
#include <pamilo/generic/benson_dual/abstract_online_vertex_enumerator.h>

namespace pamilo {

/**
 * @brief Class for online vertex enumeration
 *
 */
class OnlineVertexEnumeratorCDD : AbstractOnlineVertexEnumerator
{
public:
    OnlineVertexEnumeratorCDD() = delete;
    virtual ~OnlineVertexEnumeratorCDD();

    /**
     * @brief Indicates if an unprocessed vertex exists
     *
     */
    bool has_next() override;

    /**
     * @brief Returns next unprocessed point. Point counts as processes afterwards
     *
     * @return Point*
     */
    Point *next_vertex() override;

    /**
     * @brief Adds a new hyperplane to the vertex enumeration. A point p is on the hyperplane
     * if p*normal=rhs
     *
     * @param vertex (Is ignored)
     * @param normal Normal vector of the hyperplane
     * @param rhs Right hand side of the equation
     */
    void add_hyperplane(Point &vertex, Point &normal, double rhs) override;

    /**
     * @brief Calculates distance of a point to a hyperplane through
     * vertex*normal-rhs
     *
     * @param vertex Point
     * @param normal Normal vector of the hyperplane
     * @param rhs Right hand side of the hyperplane equation
     * @return double
     */
    double getDistance(Point &vertex, Point &normal, double rhs) override;

    /**
     * @brief Returns number of hyperplanes
     *
     * @return unsigned int
     */
    unsigned int number_of_hyperplanes()
    {
        return number_hyperplanes_;
    }

    /**
     * @brief Returns time spent in vertex enumeration (add_hyperplane method calls) since start
     *
     * @return double
     */
    double get_time()
    {
#ifdef NDEBUG
        std::cout << "End cycles: " << cycles_ << std::endl;
        std::cout << "Time: " << (cycles_ / (double)CLOCKS_PER_SEC) << std::endl;
#endif
        return cycles_ / (double)CLOCKS_PER_SEC;
    }

    /**
     * @brief Constructor
     *
     * @param initial_value
     * @param dimension
     * @param epsilon
     */
    OnlineVertexEnumeratorCDD(Point &initial_value, unsigned int dimension, double epsilon);

protected:
    int number_hyperplanes_;

    std::list<Point> unprocessed_vertices_;

    dd_MatrixPtr h_representation_;
};
}  // namespace pamilo
