include_directories(../include)
include_directories(${CPLEX_INCLUDE_DIRS})
include_directories(${CONCERT_INCLUDE_DIRS})

if(USE_CDD)
	include_directories(${CDD_INCLUDE_PATH})
endif()

set(SOURCE_FILES

	#### Header ####

	# Core
	../include/pamilo/basic/point.h
	../include/pamilo/basic/lex_point_comparator.h
	../include/pamilo/basic/abstract_solver.h

	# Benchmarks
	../include/pamilo/benchmarks/lp_parser.h

	# Benson Dual
	../include/pamilo/generic/benson_dual/abstract_online_vertex_enumerator.h
	../include/pamilo/generic/benson_dual/dual_benson_scalarizer.h
	../include/pamilo/generic/benson_dual/ove_fp_v2.h

	# Geometry Tools
	../include/pamilo/geometric/projective_geometry_utilities.h

	# Parametric Integer Linear Programming
	../include/pamilo/pilp/pilp_dual_benson.h
	../include/pamilo/pilp/ilp.h

	#### Sources ####

	# Benchmarks
	benchmarks/lp_parser.cpp

	# Benson Dual
	generic/benson_dual/ove_fp_v2.cpp
)
if(USE_CDD)
	set(SOURCE_FILES ${SOURCE_FILES}
		generic/benson_dual/ove_cdd.cpp
		../include/pamilo/generic/benson_dual/ove_cdd.h
	)
endif()

add_library(pamilo STATIC ${SOURCE_FILES})
