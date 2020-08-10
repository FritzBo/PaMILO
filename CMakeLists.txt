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
		/opt/ibm/ILOG/CPLEX_Studio128
		/opt/ibm/ILOG/CPLEX_Studio129
		/opt/ibm/ILOG/CPLEX_Studio1210
)

find_path(CPLEX_INCLUDE_DIR ilcplex/cplex.h PATHS ${CPLEX}/cplex/include)
find_path(CONCERT_INCLUDE_DIR ilconcert/iloenv.h PATHS ${CPLEX}/concert/include)
find_library(ILOCPLEX_LIB libilocplex.a PATHS ${CPLEX}/cplex/lib/x86-64_linux/static_pic)
find_library(CONCERT_LIB libconcert.a PATHS ${CPLEX}/concert/lib/x86-64_linux/static_pic)
find_library(CPLEX_LIB libcplex.a PATHS ${CPLEX}/cplex/lib/x86-64_linux/static_picu)

## CDD ##
option(USE_CDD "Whether to use cdd for vertex enumeration" OFF)
if(USE_CDD)
	find_path(CDD_INCLUDE_PATH
		NAMES cdd.h
		PATHS 	/usr/include
			/usr/local/include
			~/.local/include
			/opt/local/include
	)

	find_library(CDD_LIB
		NAMES libcdd.a
		PATHS	/usr/lib
			/usr/local/lib
			~/.local/lib
			/opt/local/lib
	)
endif()

include_directories(include)
include_directories(${CPLEX_INCLUDE_DIR})
include_directories(${CONCERT_INCLUDE_DIR})

if(USE_CDD)
	include_directories(${CDD_INCLUDE_PATH})
endif()

set(SOURCE_FILES

	#### Header ####

	# Core
	include/pamilo/basic/point.h
	include/pamilo/basic/lex_point_comparator.h
	include/pamilo/basic/abstract_solver.h

	# Benchmarks
	include/pamilo/benchmarks/lp_parser.h

	# Benson Dual
	include/pamilo/generic/benson_dual/abstract_online_vertex_enumerator.h
	include/pamilo/generic/benson_dual/dual_benson_scalarizer.h
	include/pamilo/generic/benson_dual/ove_fp_v2.h

	# Geometry Tools
	include/pamilo/geometric/projective_geometry_utilities.h

	# Parametric Integer Linear Programming
	include/pamilo/pilp/pilp_dual_benson.h
	include/pamilo/pilp/ilp.h

	#### Sources ####

	# Benchmarks
	src/benchmarks/lp_parser.cpp

	# Benson Dual
	src/generic/benson_dual/ove_fp_v2.cpp
)
if(USE_CDD)
	set(SOURCE_FILES ${SOURCE_FILES}
		src/generic/benson_dual/ove_cdd.cpp
		include/pamilo/generic/benson_dual/ove_cdd.h
	)
endif()

add_library(pamilo STATIC ${SOURCE_FILES})

include_directories(include)
include_directories(tclap)
include_directories(${CPLEX_INCLUDE_DIR})
include_directories(${CONCERT_INCLUDE_DIR})
include_directories(${PROJECT_BINARY_DIR})

if(USE_CDD)
	include_directories(${CDD_INCLUDE_PATH})
endif()

set(SOURCE_FILES
	cli/main.cpp

	#### Core
	cli/basic/modules.h
	cli/basic/modules.cpp

	#### Modules
	cli/modules/pilp_benson_module.h
	cli/modules/pilp_benson_module.cpp
)

add_executable(pamilo_cli ${SOURCE_FILES})

target_link_libraries(pamilo_cli pamilo)
target_link_libraries(pamilo_cli ${ILOCPLEX_LIB})
target_link_libraries(pamilo_cli ${CONCERT_LIB})
target_link_libraries(pamilo_cli ${CPLEX_LIB})
if(USE_CDD)
	target_link_libraries(pamilo_cli ${CDD_LIB})
endif()
if(UNIX)
	target_link_libraries(pamilo_cli pthread)
	target_link_libraries(pamilo_cli dl)
endif()

# create autogen header
configure_file("cmake/config_autogen.h.in" "${PROJECT_BINARY_DIR}/config_autogen.h")
