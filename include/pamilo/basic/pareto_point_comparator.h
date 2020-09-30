//
//  pareto_point_comparator.h
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

class ParetoDominationPointComparator {
public:
    inline ParetoDominationPointComparator(double epsilon = 0)
    : epsilon_(epsilon) { }

    inline bool operator()(const Point * p1,
                           const Point * p2) const noexcept;

    inline bool operator()(const Point& p1,
                           const Point& p2) const noexcept;

    inline static bool dominates(const Point& p1,
                                 const Point& p2,
                                 double epsilon) noexcept;

private:
    const double epsilon_;
};

inline bool ParetoDominationPointComparator::
operator()(const Point& p1,
           const Point& p2) const noexcept {

    return dominates(p1, p2, epsilon_);
}

inline bool ParetoDominationPointComparator::
operator()(const Point * p1,
           const Point * p2) const noexcept {

    return dominates(*p1, *p2, epsilon_);
}

inline bool ParetoDominationPointComparator::
dominates(const Point &p1,
          const Point &p2,
          double epsilon) noexcept {

    assert(p1.dimension() == p2.dimension());

    unsigned dimension = p1.dimension();

    bool equal = true;

    for(unsigned int i = 0; i < dimension; i++) {

        if(p1[i] - p2[i] > epsilon) {
            return false;
        }

        if(abs(p1[i] - p2[i]) > epsilon) {
            equal = false;
        }
    }

    return !equal;
}


#endif
