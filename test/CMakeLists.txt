## GTEST ##
find_path(GTEST_INCLUDE_PATH
    NAMES gtest/gtest.h
    PATHS /usr/local/
)

find_library(GTEST_MAIN
    NAMES gtest_main
    PATHS /usr/local/
)

find_library(GTEST
    NAMES gtest
    PATHS /usr/local/
)

add_subdirectory(core)
add_subdirectory(parser)
add_subdirectory(ap)
add_subdirectory(ep)
add_subdirectory(geometric)