include_directories(../../include)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})
include_directories(${GTEST_INCLUDE_PATH})

set(SOURCE_FILES
basic/lex_hungarian_test.cpp
dual_benson/ap_hungarian_benson_test.cpp
)

add_executable(ap_test ${SOURCE_FILES})

target_link_libraries(ap_test mco)
target_link_libraries(ap_test debug ${OGDF-DBG} optimized ${OGDF})
target_link_libraries(ap_test debug ${COIN-DBG} optimized ${COIN})
target_link_libraries(ap_test ${CDD})
target_link_libraries(ap_test pthread)
target_link_libraries(ap_test ${GTEST} ${GTEST_MAIN})

add_test(Assignment ap_test)