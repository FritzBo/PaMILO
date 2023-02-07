cmake_minimum_required(VERSION 2.8)

project(pamilo)

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(module_dir "${PROJECT_SOURCE_DIR}/cmake")
list(INSERT CMAKE_MODULE_PATH 0 "${module_dir}" )

set(LP_SOLVER "Gurobi" CACHE STRING "Choose the lp solver.")
set_property(CACHE LP_SOLVER PROPERTY STRINGS "Gurobi" "Cplex")

if("${LP_SOLVER}" STREQUAL "Gurobi")
    add_compile_options(-DUSE_GRB)
	set(USE_GRB ON)
elseif("${LP_SOLVER}" STREQUAL "Cplex")
    add_compile_options(-DUSE_CPLEX)
	set(USE_CPLEX ON)
endif()

enable_testing()

if(USE_CPLEX)
	set(CMAKE_CXX_FLAGS "-DIL_STD") # not needed if cplex version is >= 20
endif()

if(UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++17 -Wno-ignored-attributes")
elseif(WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
endif()

if(USE_GRB)
	find_path(
		GUROBI_INCLUDE_DIR
		NAMES gurobi_c.h
		HINTS $ENV{GUROBI_HOME}
		PATH_SUFFIXES include)
	
	find_library(
		GUROBI_LIB
		NAMES gurobi gurobi81 gurobi90 gurobi95
		HINTS $ENV{GUROBI_HOME}
		PATH_SUFFIXES lib)

	find_library(
		GUROBI_CXX_LIB
		NAMES gurobi_c++
		HINTS $ENV{GUROBI_HOME}
		PATH_SUFFIXES lib)

elseif(USE_CPLEX)

	find_path(CPLEX
		NAMES cplex
		PATHS "/usr/"
			"/usr/local/"
			"~/.local/"
			"~/opt/local"
			"/opt/ibm/ILOG/CPLEX_Studio129"
			"/opt/ibm/ILOG/CPLEX_Studio1210"
			"C:/Program Files/IBM/ILOG/CPLEX_Studio_Community1210"
			"C:/Program Files/IBM/ILOG/CPLEX_Studio1210"
			"C:/Program Files/IBM/ILOG/CPLEX_Studio129"
	)

	find_path(CPLEX_INCLUDE_DIR ilcplex/cplex.h
		PATHS "${CPLEX}/cplex/include"
	)
	find_path(CONCERT_INCLUDE_DIR ilconcert/iloenv.h
		PATHS "${CPLEX}/concert/include"
	)
	find_library(ILOCPLEX_LIB libilocplex.a
		PATHS "${CPLEX}/cplex/lib/x86-64_linux/static_pic"
			"${CPLEX}/cplex/lib/x64_windows_msvc14/stat_mda"
	)
	find_library(CONCERT_LIB libconcert.a
		PATHS "${CPLEX}/concert/lib/x86-64_linux/static_pic"
			"${CPLEX}/concert/lib/x64_windows_msvc14/stat_mda"
	)
	find_library(CPLEX_LIB libcplex.a
		PATHS "${CPLEX}/cplex/lib/x86-64_linux/static_pic"
			"${CPLEX}/cplex/lib/x64_windows_msvc14/stat_mda"
	)

endif()

## CDD ##
option(USE_CDD "Whether to use cdd for vertex enumeration" OFF)
if(USE_CDD)
	find_path(CDD_INCLUDE_PATH
		NAMES cdd.h
		PATHS 	"/usr/include/cddlib"
				"/usr/local/include/cddlib"
				"~/.local/include/cddlib"
				"/opt/local/include/cddlib"
	)

	find_library(CDD_LIB
		NAMES libcdd.so libcdd.a
		PATHS	"/usr/lib"
				"/usr/local/lib"
				"~/.local/lib"
				"/opt/local/lib"
	)
endif()

include_directories(include)
include_directories(tclap)
include_directories(cli)

if(USE_GRB)
	include_directories(${GUROBI_INCLUDE_DIR})
elseif(USE_CPLEX)
	include_directories(${CPLEX_INCLUDE_DIR})
	include_directories(${CONCERT_INCLUDE_DIR})
endif()

if(USE_CDD)
	include_directories(${CDD_INCLUDE_PATH})
endif()

set(SOURCE_FILES

	#### Header ####

	# Core
	"include/pamilo/basic/point.h"
	"include/pamilo/basic/lex_point_comparator.h"
	"include/pamilo/basic/abstract_solver.h"

	# Benchmarks
	"include/pamilo/benchmarks/lp_parser.h"

	# Benson Dual
	"include/pamilo/generic/benson_dual/abstract_online_vertex_enumerator.h"
	"include/pamilo/generic/benson_dual/dual_benson_scalarizer.h"
	"include/pamilo/generic/benson_dual/ove_fp_v2.h"

	# Geometry Tools
	"include/pamilo/geometric/projective_geometry_utilities.h"

	# Parametric Integer Linear Programming
	"include/pamilo/pilp/pilp_dual_benson.h"
	"include/pamilo/pilp/ilp.h"

	#### Sources ####

	# Benchmarks
	"src/benchmarks/lp_parser.cpp"

	# Benson Dual
	"src/generic/benson_dual/ove_fp_v2.cpp"
)

add_library(pamilo STATIC ${SOURCE_FILES})

set(SOURCE_FILES
	"cli/main.cpp"

	#### Core
	"cli/basic/modules.h"
	"cli/basic/modules.cpp"

	#### Modules
	"cli/modules/pilp_benson_module.h"
	"cli/modules/pilp_benson_module.cpp"
)

if(USE_CDD)
	set(SOURCE_FILES ${SOURCE_FILES}
		"src/generic/benson_dual/ove_cdd.cpp"
		"include/pamilo/generic/benson_dual/ove_cdd.h"
	)

	add_compile_options(-DUSE_CDD)
endif()

add_executable(pamilo_cli ${SOURCE_FILES})

target_link_libraries(pamilo_cli pamilo)

if(USE_GRB)
	target_link_libraries(pamilo_cli ${GUROBI_CXX_LIB})
	target_link_libraries(pamilo_cli ${GUROBI_LIB})
elseif(USE_CPLEX)
	target_link_libraries(pamilo_cli "${ILOCPLEX_LIB}")
	target_link_libraries(pamilo_cli "${CONCERT_LIB}")
	target_link_libraries(pamilo_cli "${CPLEX_LIB}")
endif()

if(USE_CDD)
	target_link_libraries(pamilo_cli "${CDD_LIB}")
endif()

if(UNIX)
	target_link_libraries(pamilo_cli pthread)
	target_link_libraries(pamilo_cli dl)
endif()