#pragma once
/*
 * point.h
 *
 *  Created on: 15.03.2013
 *      Author: fritz
 */

#ifndef POINT_H_
#define POINT_H_

#include <ostream>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <iostream> // FIXME

namespace mco {
    
class Point {
public:

    /// Destructor
    virtual ~Point() noexcept { delete[] values_; }
    
    inline Point() : dimension_(0), values_(nullptr) { }
    inline explicit Point(unsigned int dimension);
	inline Point(const double *values, unsigned int dimension);
    inline Point(double value, unsigned int dimension);
    inline Point(std::initializer_list<double> values);

    /// Copy constructor
	inline Point(const Point &p);
    
    /// Move constructor
    inline Point(Point&& that) noexcept;
    
    /// Copy assignment
	inline Point & operator=(const Point& that) noexcept;
    
    // Move assignment
    inline Point & operator=(Point&& that) noexcept;
    
    unsigned dimension() const noexcept { return dimension_; }

    inline Point & operator+=(const Point &that) noexcept;
    inline Point & operator-=(const Point &that) noexcept;
    
    inline Point & operator*=(double d) noexcept;

    inline Point operator-() const &;
    inline Point operator-() &&;

    inline Point operator+(Point that) const & noexcept;
    inline Point operator+(const Point& that) && noexcept;
    
    inline Point operator-(Point that) const & noexcept;
    inline Point operator-(const Point& that) && noexcept;
    
    /// Inner Product
    inline double operator*(const Point &point) const noexcept;

	double& operator[](unsigned int index) { return values_[index]; }
	const double& operator[](unsigned int index) const { return values_[index]; }

    const double* cbegin() const { return values_; }
    const double* cend() const { return values_ + dimension_; }
    
    double* begin() { return values_; }
    double* end() { return values_ + dimension_; }

	static Point * Null(const unsigned int dimension) { return new Point(0.0, dimension); }
	static Point * One(const unsigned int dimension) { return new Point(1.0, dimension); }

	friend std::ostream & operator<<(std::ostream &, const Point &);
    friend void swap(Point& p1, Point& p2);

    
private:
    unsigned int dimension_ = 0;
    double * values_ = nullptr;
};

inline Point::Point(unsigned dimension)
: dimension_(dimension) {
    values_ = dimension_ ? new double[dimension] : nullptr;
    for(unsigned int i = 0; i < dimension; ++i)
        values_[i] = 0;
}
    
inline Point::Point(const double *values, unsigned int dimension)
: dimension_(dimension) {
    values_ = dimension_ ? new double[dimension] : nullptr;
    std::copy(values, values + dimension, values_);
}
    
inline Point::Point(double value, unsigned int dimension)
: dimension_(dimension) {
    values_ = dimension_ ? new double[dimension] : nullptr;
    for(unsigned i = 0; i < dimension; ++i) {
        values_[i] = value;
    }
}
    
inline Point::Point(std::initializer_list<double> values)
: dimension_(values.size()) {
    values_ = dimension_ ? new double[values.size()] : nullptr;
    std::copy(values.begin(), values.end(), values_);
}

inline Point::Point(const Point &p) : dimension_(p.dimension_) {
    values_ = dimension_ ? new double[dimension_] : nullptr;
    std::copy(p.values_, p.values_ + dimension_, values_);
}
    
inline Point::Point(Point&& that) noexcept : Point() {
    swap(*this, that);
}
    
inline Point& Point::operator=(const Point& that) noexcept {
    
    if(that.dimension_ <= dimension_) {
        dimension_ = that.dimension_;
        std::copy(that.cbegin(), that.cend(), this->begin());
    } else {
        Point new_point(that);
        swap(*this, new_point);
    }
    
    return *this;
}
    
inline Point& Point::operator=(Point&& that) noexcept {
    swap(*this, that);
    return *this;
}

    
inline Point & Point::operator+=(const Point &p) noexcept {
    assert(dimension_ == p.dimension_);
    
    for(unsigned int i = 0; i < dimension_; ++i)
        values_[i] += p[i];
        
    return *this;
}

inline Point & Point::operator-=(const Point &p) noexcept {
    assert(dimension_ == p.dimension_);
    
    for(unsigned int i = 0; i < dimension_; ++i)
        values_[i] -= p[i];
        
    return *this;
}
    
inline Point & Point::operator*=(double d) noexcept {
    for(unsigned int i = 0; i < dimension_; ++i)
        values_[i] *= d;
        
    return *this;
}
    
inline Point Point::operator-() const & {
    Point result(*this);
    result *= -1;
    return result;
}
    
inline Point Point::operator-() && {
    *this *= -1;
    return std::move(*this);
}
    
inline Point Point::operator+(Point that) const & noexcept {
    assert(dimension_ == that.dimension_);
    return std::move(that += *this);
}

inline Point Point::operator+(const Point& that) && noexcept {
    assert(dimension_ == that.dimension_);
    return std::move(*this += that);
}
    
inline Point Point::operator-(Point that) const & noexcept {
    assert(dimension_ == that.dimension_);
    
    for(unsigned i = 0; i < dimension_; ++i) {
        that[i] = this->operator[](i) - that[i];
    }
    
    return std::move(that);
}
    
inline Point Point::operator-(const Point& that) && noexcept {
    assert(dimension_ == that.dimension_);
    return std::move(*this -= that);
}
    
inline double Point::operator*(const Point &point) const noexcept {
    assert(dimension_ == point.dimension_);
    
    double sum = 0;
    unsigned i;
    for(i = 1; i < dimension_; i += 2) {
        sum += (*this)[i - 1] * point[i - 1];
        sum += (*this)[i] * point[i];
    }
    
    return i == dimension_ + 1 ? sum : sum + (*this)[dimension_ - 1] * point[dimension_-1];
}

inline Point operator*(Point p, double d) {
    return std::move(p *= d);
}
    
inline Point operator*(double d, Point p) {
    return std::move(p *= d);
}
    
inline std::ostream & operator<<(std::ostream &os, const Point &point) {
    os << "(";
    
    if(point.dimension_ != 0) {
        for(unsigned int i = 0; i < point.dimension_ - 1; ++i) {
            os << point[i] << ", ";
        }
        
        os << point[point.dimension_ - 1];
    }
    
    os << ")";
    
    return os;
}
    
inline void swap(Point& p1, Point& p2) {
    using std::swap;
    
    swap(p1.values_,    p2.values_);
    swap(p1.dimension_, p2.dimension_);
}

#include "lex_point_comparator.h"
#include "equality_point_comparator.h"
#include "componentwise_point_comparator.h"
#include "pareto_point_comparator.h"
    
} // namespace mco


#endif /* POINT_H_ */
