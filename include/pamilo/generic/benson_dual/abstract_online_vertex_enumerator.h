/**
 * @file abstract_vertex_enumerator.h
 * @author Fritz Bökler
 * @brief
 * @date 08.12.2013
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#pragma once

#include <ctime>
#include <list>
#include <map>
#include <queue>
#include <vector>

#include <pamilo/basic/point.h>

namespace pamilo {

/**
 * @brief Interface for the Online Vertex Enumerator
 *
 */
class AbstractOnlineVertexEnumerator
{
public:
    AbstractOnlineVertexEnumerator() = delete;
    virtual ~AbstractOnlineVertexEnumerator() = default;

    /**
     * @brief Returns the cpu time
     *
     * @return double
     */
    double get_time()
    {
        return cycles_ / (double)CLOCKS_PER_SEC;
    }

    /**
     * @brief Construct a new Abstract Online Vertex Enumerator object
     *
     * @param dimension
     * @param epsilon
     */
    AbstractOnlineVertexEnumerator(unsigned int dimension, double epsilon)
        : dimension_(dimension)
        , epsilon_(epsilon)
        , cycles_(0)
    {
    }

    /**
     * @brief Indicates whether an unprocessed vertex exists
     *
     */
    virtual inline bool has_next() = 0;

    /**
     * @brief Returns next unprocessed point. Point counts as processes afterwards
     *
     * @return Point*
     */
    virtual inline Point *next_vertex() = 0;

    /**
     * @brief Adds a new hyperplane to the vertex enumeration. A point p is on the hyperplane
     * if p*normal=rhs
     *
     * @param vertex
     * @param normal Normal vector of the hyperplane
     * @param rhs Right hand side of the equation
     */
    virtual void add_hyperplane(Point &vertex, Point &normal, double rhs) = 0;

    /**
     * @brief Returns number of hyperplanes
     *
     * @return unsigned int
     */
    virtual unsigned int number_of_hyperplanes() = 0;

    /**
     * @brief Calculates distance of a point to a hyperplane through
     * vertex*normal-rhs
     *
     * @param vertex Point
     * @param normal Normal vector of the hyperplane
     * @param rhs Right hand side of the hyperplane equation
     * @return double
     */
    virtual double getDistance(Point &vertex, Point &normal, double rhs) = 0;

protected:
    /**
     * @brief dimension of the space
     *
     */
    const unsigned int dimension_;

    /**
     * @brief epsilon for floating point comparisons
     *
     */
    const double epsilon_;

    /**
     * @briefCPU time spent on vertex enumeration
     *
     */
    clock_t cycles_;
};
}  // namespace pamilo
