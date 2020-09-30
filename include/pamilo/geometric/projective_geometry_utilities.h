//
//  projective_geometry_utilities.h
//  pamilo
//
//  Created by Fritz Bökler on 11.04.14.
//
//  This file is distributed for academics only
//  under the terms of an MIT license based license,
//  a copy of which can be found in the file LICENSE-academic.txt.
//
//

#pragma once

#include <type_traits>

#include <pamilo/basic/point.h>

namespace pamilo {

class ProjectiveGeometry {
public:
    template<typename PointType,
    typename = typename std::enable_if<std::is_convertible<PointType, Point>::value>::type>
    inline static PointType to_projective(const Point& point);

    template<typename PointType,
    typename = typename std::enable_if<std::is_convertible<PointType, Point>::value>::type>
    inline static PointType from_projective(const Point& Point);

    inline static void normalize_projective(Point& point);
};

template<typename PointType, typename>
inline PointType ProjectiveGeometry::
to_projective(const Point &point) {

    unsigned dimension = point.dimension();
    PointType new_proj_point(dimension + 1);

	for(unsigned int i = 0; i < dimension; ++i)
		new_proj_point[i] = point[i];

	new_proj_point[dimension] = 1;

	return new_proj_point;
}

template<typename PointType, typename>
inline PointType ProjectiveGeometry::
from_projective(const Point& point) {

    unsigned projective_dimension = point.dimension();
    PointType new_point(projective_dimension - 1);

    double normalization = point[projective_dimension - 1];

    for(unsigned i = 0; i < projective_dimension - 1; ++i) {
        new_point[i] = point[i] / normalization;
    }

    return new_point;
}


inline void ProjectiveGeometry::
normalize_projective(Point& projective_point) {

    unsigned dimension = projective_point.dimension();

    for(unsigned int i = 0; i < dimension; ++i)
        projective_point[i] = projective_point[i] / projective_point[dimension - 1];

    projective_point[dimension - 1] = 1;
}
}

