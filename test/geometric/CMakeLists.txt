include_directories(../../include)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})
include_directories(${GTEST_INCLUDE_PATH})
include_directories(${CDD_INCLUDE_PATH})

set(SOURCE_FILES
ove_fp_v2_test.cpp
)

add_executable(geometry_test ${SOURCE_FILES})

target_link_libraries(geometry_test mco)
target_link_libraries(geometry_test debug ${OGDF-DBG} optimized ${OGDF})
target_link_libraries(geometry_test debug ${COIN-DBG} optimized ${COIN})
target_link_libraries(geometry_test ${CDD})
target_link_libraries(geometry_test pthread)
target_link_libraries(geometry_test ${GTEST} ${GTEST_MAIN})

add_test(Geometry_Tests geometry_test)
