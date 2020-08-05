cmake_minimum_required(VERSION 2.6)

project(pamilo)

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

find_path(CPLEX_INCLUDE_DIRS ilcplex/cplex.h PATHS ${CPLEX}/cplex/include)
find_path(CONCERT_INCLUDE_DIRS ilconcert/iloenv.h PATHS ${CPLEX}/concert/include)
find_library(ILOCPLEX_LIB libilocplex.a PATHS ${CPLEX}/cplex/lib/x86-64_linux/static_pic)
find_library(CONCERT_LIB libconcert.a PATHS ${CPLEX}/concert/lib/x86-64_linux/static_pic)
find_library(CPLEX_LIB libcplex.a PATHS ${CPLEX}/cplex/lib/x86-64_linux/static_pic)

## CDD ##
option(USE_CDD "Whether to use cdd for vertex enumeration" OFF)
if(USE_CDD)
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
endif()

add_subdirectory(pamilo)
add_subdirectory(cli)

# create autogen header
configure_file("cmake/config_autogen.h.in" "${PROJECT_BINARY_DIR}/config_autogen.h")
