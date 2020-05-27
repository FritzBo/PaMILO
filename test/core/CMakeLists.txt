include_directories(../../include)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})
include_directories(${GTEST_INCLUDE_PATH})

set(SOURCE_FILES
Point_test.cpp
)

add_executable(core_test ${SOURCE_FILES})

target_link_libraries(core_test mco)
target_link_libraries(core_test debug ${OGDF-DBG} optimized ${OGDF})
target_link_libraries(core_test debug ${COIN-DBG} optimized ${COIN})
target_link_libraries(core_test pthread)
target_link_libraries(core_test ${GTEST} ${GTEST_MAIN})

add_test(Core core_test)
