cmake_minimum_required(VERSION 2.6)

project(mco)

set(module_dir "${PROJECT_SOURCE_DIR}/cmake")
list(INSERT CMAKE_MODULE_PATH 0 "${module_dir}" )

enable_testing()

set(CMAKE_CXX_FLAGS
--std=c++17)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DOGDF_DEBUG")

## OGDF and COIN ##
set(OGDF_CONFIG_SEARCH_PATH CACHE PATH "")

find_package(OGDF CONFIG PATHS ${OGDF_CONFIG_SEARCH_PATH})

find_library(OSI
	NAMES OSI
		libOsiCpx.so
    PATHS /usr/local/
	~/.local/
)

## CPLEX
find_library(CPLEX_LIB_PATH
	NAMES cplex)

find_path(OSI_INCLUDE_DIRS
    NAMES coin
    PATHS /usr/local/
	~/.local/
)


## CDD ##
find_path(CDD_INCLUDE_PATH
	NAMES cdd.h
	PATHS 	/usr/
		/usr/local/
	~/.local/
		/opt/local/
)

find_library(CDD
	NAMES libcdd.a
	PATHS	/usr/
		/usr/local
	~/.local/lib
		/opt/local
)

add_subdirectory(mco)
add_subdirectory(cli)
