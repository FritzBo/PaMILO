//
//  ove_fp_v2_test.cpp
//  mco
//
//  Created by Fritz Bökler on 13.04.14.
//
//

#include <vector>
#include <list>

using std::vector;
using std::list;

#include <gtest/gtest.h>

using ::testing::Combine;
using ::testing::Values;

#include <mco/basic/point.h>
#include <mco/generic/benson_dual/ove_fp_v2.h>

using mco::Point;
using mco::GraphlessOVE;

TEST(GraphlessOVETest, UnitCubeTest) {
    std::list<Point> unit_cube_initial_inequalities = {
        Point({1, 0, 0, -1}),
        Point({0, 1, 0, -1}),
        Point({0, 0, 1, -1}),
    };
    
    std::list<Point> unit_cube_inequalities = {
        Point({-1, 0, 0}),
        Point({0, -1, 0}),
        Point({0, 0, -1}),
    };
    
    std::list<Point> unit_cube_initial_extreme_points = {
        Point({-1, -1, -1})
    };
    
    list<Point> empty_range = {};
    
    GraphlessOVE ove(3,
                     unit_cube_initial_extreme_points.cbegin(),
                     unit_cube_initial_extreme_points.cend(),
                     empty_range.cbegin(),
                     empty_range.cend(),
                     unit_cube_initial_inequalities.cbegin(),
                     unit_cube_initial_inequalities.cend());
    
    while(ove.has_next()) {
        Point* candidate = ove.next_vertex();
        
        if(!unit_cube_inequalities.empty()) {
            ove.add_hyperplane(*candidate,
                               unit_cube_inequalities.back(),
                               -1);
            unit_cube_inequalities.pop_back();
        }
        
//        cout << *candidate << endl;
        
        delete candidate;
    }
}