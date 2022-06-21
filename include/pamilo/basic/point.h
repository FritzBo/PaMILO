/**
 * @file point.h
 * @author Fritz Bökler
 * @brief
 * @date 15.03.2013
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */
#pragma once

#include <cassert>
#include <cmath>
#include <iostream>  // FIXME
#include <ostream>
#include <stdexcept>

namespace pamilo {

/**
 * @brief Represents a point of dimension dimension_
 *
 */
class Point
{
public:
    /**
     * @brief Destroy the Point object
     *
     */
    virtual ~Point() noexcept
    {
        delete[] values_;
    }

    /**
     * @brief Constructs a point of dimension 0
     *
     */
    inline Point()
        : dimension_(0)
        , values_(nullptr)
    {
    }

    /**
     * @brief Construct a new Point object. All values are initially 0
     *
     * @param dimension Number of dimensions
     */
    inline explicit Point(unsigned int dimension);

    /**
     * @brief Construct a new Point object with values
     *
     * @param values Pointer to an array of values. Copies the first dimension entries of values
     * into the point
     * @param dimension Number of dimensions
     */
    inline Point(const double *values, unsigned int dimension);

    /**
     * @brief Construct a new Point object and initialises values to value
     *
     * @param value Initial value for all values
     * @param dimension Number of dimensions
     */
    inline Point(double value, unsigned int dimension);

    /**
     * @brief Construct a new Point object from a list of values. Dimension is size of the list
     *
     * @param values Initial values
     */
    inline Point(const std::initializer_list<double> values);

    /**
     * @brief Copy constructor
     *
     * @param p
     */
    inline Point(const Point &p);

    /**
     * @brief Move constructor
     *
     * @param that
     */
    inline Point(Point &&that) noexcept;

    /**
     * @brief Copy assignment
     *
     * @param that
     * @return Point&
     */
    inline Point &operator=(const Point &that) noexcept;

    /**
     * @brief Move assignment
     *
     * @param that
     * @return Point&
     */
    inline Point &operator=(Point &&that) noexcept;

    /**
     * @brief dimension getter
     *
     * @return unsigned
     */
    unsigned dimension() const noexcept
    {
        return dimension_;
    }

    /**
     * @brief Adds another point componentwise to this point
     *
     * @param that
     * @return Point&
     */
    inline Point &operator+=(const Point &that) noexcept;

    /**
     * @brief Substracts another point componentwise from this point
     *
     * @param that
     * @return Point&
     */
    inline Point &operator-=(const Point &that) noexcept;

    /**
     * @brief Scalar multiplication on this point
     *
     * @param d
     * @return Point&
     */
    inline Point &operator*=(double d) noexcept;

    /**
     * @brief Point substraction
     *
     * @return Point
     */
    inline Point operator-() const &;

    /**
     * @brief Point substraction
     *
     * @return Point
     */
    inline Point operator-() &&;

    /**
     * @brief Point addition
     *
     * @param that
     * @return Point
     */
    inline Point operator+(Point that) const &noexcept;

    /**
     * @brief Point addition
     *
     * @param that
     * @return Point
     */
    inline Point operator+(const Point &that) &&noexcept;

    /**
     * @brief Point substraction
     *
     * @param that
     * @return Point
     */
    inline Point operator-(Point that) const &noexcept;

    /**
     * @brief Point substraction
     *
     * @param that
     * @return Point
     */
    inline Point operator-(const Point &that) &&noexcept;

    /**
     * @brief Calculates inner product
     *
     * @param point
     * @return double
     */
    inline double operator*(const Point &point) const noexcept;

    /**
     * @brief Access to values
     *
     * @param index
     * @return double&
     */
    double &operator[](unsigned int index)
    {
        return values_[index];
    }

    /**
     * @brief Const access to values
     *
     * @param index
     * @return const double&
     */
    const double &operator[](unsigned int index) const
    {
        return values_[index];
    }

    /**
     * @brief Const pointer to begin of values array
     *
     * @return const double*
     */
    const double *cbegin() const
    {
        return values_;
    }

    /**
     * @brief Const pointer to end of values array
     *
     * @return const double*
     */
    const double *cend() const
    {
        return values_ + dimension_;
    }

    /**
     * @brief Pointer to begin of values array
     *
     * @return double*
     */
    double *begin()
    {
        return values_;
    }

    /**
     * @brief Pointer to end of values array
     *
     * @return double*
     */
    double *end()
    {
        return values_ + dimension_;
    }

    /**
     * @brief Constructs new point of dimension dimensions and initial values 0 on the heap
     *
     * @param dimension
     * @return Point*
     */
    static Point *Null(const unsigned int dimension)
    {
        return new Point(0.0, dimension);
    }

    /**
     * @brief Constructs new point of dimension dimensions and initial values 1 on the heap
     *
     * @param dimension
     * @return Point*
     */
    static Point *One(const unsigned int dimension)
    {
        return new Point(1.0, dimension);
    }

    /**
     * @brief Prints the entries of values in order on the ostream
     *
     * @return std::ostream& same as input ostream&)
     */
    friend std::ostream &operator<<(std::ostream &, const Point &);

    /**
     * @brief Swap implementation
     *
     * @param p1
     * @param p2
     */
    friend void swap(Point &p1, Point &p2);

private:
    /**
     * @brief Dimension of the point
     *
     */
    unsigned int dimension_ = 0;
    /**
     * @brief Values of each dimension
     *
     */
    double *values_ = nullptr;
};

inline Point::Point(unsigned dimension)
    : dimension_(dimension)
{
    values_ = dimension_ ? new double[dimension] : nullptr;
    for (unsigned int i = 0; i < dimension; ++i)
        values_[i] = 0;
}

inline Point::Point(const double *values, unsigned int dimension)
    : dimension_(dimension)
{
    values_ = dimension_ ? new double[dimension] : nullptr;
    std::copy(values, values + dimension, values_);
}

inline Point::Point(double value, unsigned int dimension)
    : dimension_(dimension)
{
    values_ = dimension_ ? new double[dimension] : nullptr;
    for (unsigned i = 0; i < dimension; ++i)
    {
        values_[i] = value;
    }
}

inline Point::Point(const std::initializer_list<double> values)
    : dimension_(values.size())
{
    values_ = dimension_ ? new double[values.size()] : nullptr;
    std::copy(values.begin(), values.end(), values_);
}

inline Point::Point(const Point &p)
    : dimension_(p.dimension_)
{
    values_ = dimension_ ? new double[dimension_] : nullptr;
    std::copy(p.values_, p.values_ + dimension_, values_);
}

inline Point::Point(Point &&that) noexcept
    : Point()
{
    swap(*this, that);
}

inline Point &Point::operator=(const Point &that) noexcept
{
    if (that.dimension_ <= dimension_)
    {
        dimension_ = that.dimension_;
        std::copy(that.cbegin(), that.cend(), this->begin());
    }
    else
    {
        Point new_point(that);
        swap(*this, new_point);
    }

    return *this;
}

inline Point &Point::operator=(Point &&that) noexcept
{
    swap(*this, that);
    return *this;
}

inline Point &Point::operator+=(const Point &p) noexcept
{
    assert(dimension_ == p.dimension_);

    for (unsigned int i = 0; i < dimension_; ++i)
        values_[i] += p[i];

    return *this;
}

inline Point &Point::operator-=(const Point &p) noexcept
{
    assert(dimension_ == p.dimension_);

    for (unsigned int i = 0; i < dimension_; ++i)
        values_[i] -= p[i];

    return *this;
}

inline Point &Point::operator*=(double d) noexcept
{
    for (unsigned int i = 0; i < dimension_; ++i)
        values_[i] *= d;

    return *this;
}

inline Point Point::operator-() const &
{
    Point result(*this);
    result *= -1;
    return result;
}

inline Point Point::operator-() &&
{
    *this *= -1;
    return std::move(*this);
}

inline Point Point::operator+(Point that) const &noexcept
{
    assert(dimension_ == that.dimension_);
    return std::move(that += *this);
}

inline Point Point::operator+(const Point &that) &&noexcept
{
    assert(dimension_ == that.dimension_);
    return std::move(*this += that);
}

inline Point Point::operator-(Point that) const &noexcept
{
    assert(dimension_ == that.dimension_);

    for (unsigned i = 0; i < dimension_; ++i)
    {
        that[i] = this->operator[](i) - that[i];
    }

    return std::move(that);
}

inline Point Point::operator-(const Point &that) &&noexcept
{
    assert(dimension_ == that.dimension_);
    return std::move(*this -= that);
}

inline double Point::operator*(const Point &point) const noexcept
{
    assert(dimension_ == point.dimension_);

    double sum = 0;
    unsigned i;
    for (i = 1; i < dimension_; i += 2)
    {
        sum += (*this)[i - 1] * point[i - 1];
        sum += (*this)[i] * point[i];
    }

    return i == dimension_ + 1 ? sum : sum + (*this)[dimension_ - 1] * point[dimension_ - 1];
}

inline Point operator*(Point p, double d)
{
    return std::move(p *= d);
}

inline Point operator*(double d, Point p)
{
    return std::move(p *= d);
}

inline std::ostream &operator<<(std::ostream &os, const Point &point)
{
    os << "";

    if (point.dimension_ != 0)
    {
        for (unsigned int i = 0; i < point.dimension_ - 1; ++i)
        {
            os << point[i] << " ";
        }

        os << point[point.dimension_ - 1];
    }

    os << "";

    return os;
}

inline void swap(Point &p1, Point &p2)
{
    using std::swap;

    swap(p1.values_, p2.values_);
    swap(p1.dimension_, p2.dimension_);
}

#include "componentwise_point_comparator.h"
#include "equality_point_comparator.h"
#include "lex_point_comparator.h"
#include "pareto_point_comparator.h"
}  // namespace pamilo
