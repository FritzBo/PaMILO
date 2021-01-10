//
//  equality_point_comparator.h
//  pamilo
//
//  Created by Fritz Bökler on 31.03.14.
//
//  This file is distributed for academics only
//  under the terms of an MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
//
//

#pragma once

class EqualityPointComparator {
public:
    inline EqualityPointComparator(double epsilon)
    : epsilon_(epsilon) { }

    inline bool operator()(const Point * point1,
                           const Point * point2) const noexcept;

    inline bool operator()(const Point& point1,
                           const Point& point2) const noexcept;

    inline static bool is_equal(const Point& point1,
                                const Point& point2,
                                double epsilon) noexcept;

private:
    const double epsilon_;
};

inline bool EqualityPointComparator::
operator()(const Point * p1,
           const Point * p2) const noexcept {

    return is_equal(*p1, *p2, epsilon_);
}

inline bool EqualityPointComparator::
operator()(const Point& p1,
           const Point& p2) const noexcept {

    return is_equal(p1, p2, epsilon_);
}


inline bool EqualityPointComparator::
is_equal(const Point& p1,
         const Point& p2,
         double epsilon) noexcept {

    if(p1.dimension() != p2.dimension()) {
        return false;
    }

    unsigned dimension = p1.dimension();

    for(unsigned int i = 0; i < dimension; ++i) {
		if(std::abs(p1[i] - p2[i]) > epsilon) {
			return false;
        }
    }

	return true;
}

