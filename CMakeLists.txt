cmake_minimum_required(VERSION 2.6)

project(mco)

set(module_dir "${PROJECT_SOURCE_DIR}/cmake")
list(INSERT CMAKE_MODULE_PATH 0 "${module_dir}" )

enable_testing()

set(CMAKE_CXX_FLAGS
	"--std=c++17 -DIL_STD -Wno-ignored-attributes")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")

## CPLEX
find_path(CPLEX
	NAMES cplex
	PATHS /usr/
		/usr/local/
		~/.local/
		~/opt/local
)

set(CPLEX_INCLUDE_DIRS ${CPLEX}/cplex/include)
set(CONCERT_INCLUDE_DIRS ${CPLEX}/concert/include)
set(ILOCPLEX_LIB_PATH ${CPLEX}/cplex/lib/x86-64_linux/static_pic/libilocplex.a)
set(CONCERT_LIB_PATH ${CPLEX}/concert/lib/x86-64_linux/static_pic/libconcert.a)
set(CPLEX_LIB_PATH ${CPLEX}/cplex/lib/x86-64_linux/static_pic/libcplex.a)

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
