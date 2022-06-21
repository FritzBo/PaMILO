/**
 * @file lex_point_comparator.h
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
 * @brief Compares two points. Returns ture if first point is lexicographic smaller. With strict
 * also returns true if points are equal
 *
 */
class LexPointComparator
{
public:
    inline LexPointComparator(double epsilon, bool strict = true)
        : epsilon_(epsilon)
        , strict_(strict)
    {
    }

    inline bool operator()(const Point *p1, const Point *p2) const noexcept;

    inline bool operator()(const Point &p1, const Point &p2) const noexcept;

    inline static bool is_lex_le(const Point &p1, const Point &p2, double epsilon) noexcept;

    inline static bool is_lex_leq(const Point &p1, const Point &p2, double epsilon) noexcept;

private:
    const double epsilon_;
    const bool strict_;
};

inline bool LexPointComparator::operator()(const Point &p1, const Point &p2) const noexcept
{
    if (strict_)
    {
        return is_lex_le(p1, p2, epsilon_);
    }
    else
    {
        return is_lex_leq(p1, p2, epsilon_);
    }
}

inline bool LexPointComparator::operator()(const Point *point1, const Point *point2) const noexcept
{
    return operator()(*point1, *point2);
}

inline bool LexPointComparator::is_lex_le(const Point &p1, const Point &p2, double epsilon) noexcept
{
    assert(p1.dimension() == p2.dimension());

    unsigned dimension = p1.dimension();

    for (unsigned int i = 0; i < dimension; ++i)
    {
        if (p1[i] - p2[i] < -epsilon)
        {
            return true;
        }
        else if (p1[i] - p2[i] > epsilon)
        {
            return false;
        }
    }
    return false;
}

inline bool LexPointComparator::is_lex_leq(const Point &p1, const Point &p2,
                                           double epsilon) noexcept
{
    assert(p1.dimension() == p2.dimension());

    unsigned dimension = p1.dimension();

    for (unsigned int i = 0; i < dimension; ++i)
    {
        if (p1[i] - p2[i] < -epsilon)
        {
            return true;
        }
        else if (p1[i] - p2[i] > epsilon)
        {
            return false;
        }
    }
    return true;
}

class LexWeightedPointComparator
{
public:
    inline LexWeightedPointComparator(Point weighting, double epsilon, bool strict = true)
        : epsilon_(epsilon)
        , strict_(strict)
        , weighting_(weighting)
    {
    }

    inline bool operator()(const Point *p1, const Point *p2) const noexcept;

    inline bool operator()(const Point &p1, const Point &p2) const noexcept;

private:
    const double epsilon_;
    const bool strict_;

    Point weighting_;
};

inline bool LexWeightedPointComparator::operator()(const Point &p1, const Point &p2) const noexcept
{
    double weight1 = p1 * weighting_;
    double weight2 = p2 * weighting_;

    if (weight1 < weight2)
    {
        return true;
    }
    else if (weight2 < weight1)
    {
        return false;
    }

    if (strict_)
    {
        return LexPointComparator::is_lex_le(p1, p2, epsilon_);
    }
    else
    {
        return LexPointComparator::is_lex_leq(p1, p2, epsilon_);
    }
}

inline bool LexWeightedPointComparator::operator()(const Point *p1, const Point *p2) const noexcept
{
    return operator()(*p1, *p2);
}
