//
//  componentwise_point_comparator.h
//  pamilo
//
//  Created by Fritz Bökler on 31.03.14.
//
//  This file is distributed under the terms of
//
//  the GNU General Public License v3,
//  a copy of which can be found in the file LICENCE-GPLv3.txt
//
//  OR
//
//  for academics, a MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
//
//

#pragma once

class ComponentwisePointComparator {
public:
    inline ComponentwisePointComparator(double epsilon = 0, bool strict = true)
    : epsilon_(epsilon), strict_(strict) { }

    inline bool operator()(const Point * point1,
                           const Point * point2) const noexcept;

    inline bool operator()(const Point& p1,
                           const Point& p2) const noexcept;

    inline static bool is_le(const Point& p1,
                             const Point& p2,
                             double epsilon) noexcept;

    inline static bool is_leq(const Point& p1,
                              const Point& p2,
                              double epsilon) noexcept;

private:
    const double epsilon_;
    const bool strict_;
};

inline bool ComponentwisePointComparator::
operator()(const Point* p1,
           const Point* p2) const noexcept {

    if(strict_) {
        return is_le(*p1, *p2, epsilon_);
    } else {
        return is_leq(*p1, *p2, epsilon_);
    }
}

inline bool ComponentwisePointComparator::
operator()(const Point& p1,
           const Point& p2) const noexcept {

    if(strict_) {
        return is_le(p1, p2, epsilon_);
    } else {
        return is_leq(p1, p2, epsilon_);
    }
}

inline bool ComponentwisePointComparator::
is_le(const Point &p1,
      const Point &p2,
      double epsilon) noexcept {

    assert(p1.dimension() == p2.dimension());

    unsigned dimension = p1.dimension();

    for(unsigned int i = 0; i < dimension; ++i) {
        if(p1[i] - p2[i] >= - epsilon) {
            return false;
        }
    }

    return true;
}


inline bool ComponentwisePointComparator::
is_leq(const Point &p1,
       const Point &p2,
       double epsilon) noexcept {

    assert(p1.dimension() == p2.dimension());

    unsigned dimension = p1.dimension();

    for(unsigned int i = 0; i < dimension; ++i) {
        if(p1[i] - p2[i] > epsilon) {
            return false;
        }
    }

    return true;
}

#endif
