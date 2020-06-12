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

#find_package(LIBCOIN)
#find_library(OGDF
#    NAMES libOGDF.a
#    PATHS /usr/local/
#	~/.local/
#)
#
#find_library(OGDF-DBG
#	NAMES OGDF-DBG
#		libOGDF.a
#    PATHS /usr/local/
#		~/.local/
#)
#
#find_path(OGDF_INCLUDE_PATH
#    NAMES ogdf/basic/Graph.h
#    PATHS /usr/local/
#	~/.local/
#)

find_library(OSI
	NAMES OSI
		libOsiCpx.so
    PATHS /usr/local/
	~/.local/
)
#
#find_library(COIN-DBG
#    NAMES COIN-dbg
#    PATHS /usr/local/
#	~/.local/
#)
#

## CPLEX
find_library(CPLEX_LIB_PATH
	NAMES cplex)

find_path(OSI_INCLUDE_DIRS
    NAMES coin
    PATHS /usr/local/
	~/.local/
)

### GUROBI ##
#
#find_path(GUROBI_INCLUDE_PATH
#    NAMES gurobi_c++.h
#    PATHS /Library/*/mac64/include/
#          /usr/local/
#		~/.local/
#		/home/wirkom/Downloads/gurobi902/linux64/include
#)
#
#find_library(GUROBI1
#    NAMES gurobi_c++
#		libgurobi_c++.a
#    PATHS /Library/*/mac64/lib/
#          /usr/local/
#		~/.local/
#		/home/wirkom/Downloads/gurobi902/linux64/lib
#)
#          
#find_library(GUROBI2
#    NAMES gurobi56
#    PATHS /Library/*/mac64/lib/
#          /usr/local/
#	~/.local/
#)

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
#add_subdirectory(test)
add_subdirectory(cli)
