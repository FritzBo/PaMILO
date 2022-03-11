/**
 * @file projective_geometry_utilities.h
 * @author Fritz Bökler
 * @brief
 * @date 11.04.2014
 *
 * This file is distributed for academics only under the terms of an MIT license based license, a
 * copy of which can be found in the file LICENSE-academic.txt.
 *
 */

#pragma once

#include <type_traits>

#include <pamilo/basic/point.h>

namespace pamilo {

/**
 * @brief Provides static utility function to project points
 *
 */
class ProjectiveGeometry
{
public:
    /**
     * @brief Converts a point to a new projective point. Dimension of the new point is dimension of
     * old point + 1. First entries are copied, last is set to 1.
     *
     * @tparam PointType, point must be convertible into it
     * @param point Point to convert
     * @return PointType
     */
    template <typename PointType, typename = typename std::enable_if<
                                      std::is_convertible<PointType, Point>::value>::type>
    inline static PointType to_projective(const Point &point);

    /**
     * @brief Converts a projective point to a new non-projective point. Dimension of the new point
     * is dimension of old point - 1. First entries are copied and each divided by last entry.
     *
     * @tparam PointType, point must be convertible into it
     * @param point Point to convert
     * @return PointType
     */
    template <typename PointType, typename = typename std::enable_if<
                                      std::is_convertible<PointType, Point>::value>::type>
    inline static PointType from_projective(const Point &Point);

    /**
     * @brief Normalizes a projective point by dividing each entry by the last entry
     *
     * @param point Point to normalize
     */
    inline static void normalize_projective(Point &point);
};

template <typename PointType, typename>
inline PointType ProjectiveGeometry::to_projective(const Point &point)
{
    unsigned dimension = point.dimension();
    PointType new_proj_point(dimension + 1);

    for (unsigned int i = 0; i < dimension; ++i)
        new_proj_point[i] = point[i];

    new_proj_point[dimension] = 1;

    return new_proj_point;
}

template <typename PointType, typename>
inline PointType ProjectiveGeometry::from_projective(const Point &point)
{
    unsigned projective_dimension = point.dimension();
    PointType new_point(projective_dimension - 1);

    double normalization = point[projective_dimension - 1];

    for (unsigned i = 0; i < projective_dimension - 1; ++i)
    {
        new_point[i] = point[i] / normalization;
    }

    return new_point;
}

inline void ProjectiveGeometry::normalize_projective(Point &projective_point)
{
    unsigned dimension = projective_point.dimension();

    for (unsigned int i = 0; i < dimension; ++i)
        projective_point[i] = projective_point[i] / projective_point[dimension - 1];

    projective_point[dimension - 1] = 1;
}
}  // namespace pamilo
