/**
 * @file pareto_point_comparator.h
 * @author Fritz Bökler
 * @brief
 * @date 31.3.2014
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#pragma once

/**
 * @brief Compares two points. Returns true if the first point dominates the second
 *
 */
class ParetoDominationPointComparator
{
public:
    inline ParetoDominationPointComparator(double epsilon)
        : epsilon_(epsilon)
    {
    }

    inline bool operator()(const Point *p1, const Point *p2) const noexcept;

    inline bool operator()(const Point &p1, const Point &p2) const noexcept;

    inline static bool dominates(const Point &p1, const Point &p2, double epsilon) noexcept;

private:
    const double epsilon_;
};

inline bool ParetoDominationPointComparator::operator()(const Point &p1,
                                                        const Point &p2) const noexcept
{
    return dominates(p1, p2, epsilon_);
}

inline bool ParetoDominationPointComparator::operator()(const Point *p1,
                                                        const Point *p2) const noexcept
{
    return dominates(*p1, *p2, epsilon_);
}

inline bool ParetoDominationPointComparator::dominates(const Point &p1, const Point &p2,
                                                       double epsilon) noexcept
{
    assert(p1.dimension() == p2.dimension());

    unsigned dimension = p1.dimension();

    bool equal = true;

    for (unsigned int i = 0; i < dimension; i++)
    {
        if (p1[i] - p2[i] > epsilon)
        {
            return false;
        }

        if (abs(p1[i] - p2[i]) > epsilon)
        {
            equal = false;
        }
    }
    return !equal;
}
