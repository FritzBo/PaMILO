include_directories(../include)
include_directories(${OSI_INCLUDE_DIRS})
include_directories(${OGDF_INCLUDE_DIRS})

set(SOURCE_FILES

#### Header ####

# Core
../include/mco/basic/point.h
#../include/mco/basic/pareto_point_comparator.h
#../include/mco/basic/componentwise_point_comparator.h
#../include/mco/basic/equality_point_comparator.h
../include/mco/basic/lex_point_comparator.h
../include/mco/basic/abstract_solver.h
#../include/mco/basic/weight_function_adaptors.h
#../include/mco/basic/utility.h

# Benchmarks
../include/mco/benchmarks/lp_parser.h

# Benson Dual
../include/mco/generic/benson_dual/abstract_online_vertex_enumerator.h
../include/mco/generic/benson_dual/dual_benson_scalarizer.h
../include/mco/generic/benson_dual/ove_cdd.h
../include/mco/generic/benson_dual/ove_node_lists.h
../include/mco/generic/benson_dual/ove_edge_lists.h
../include/mco/generic/benson_dual/ove_fp_v2.h

# Geometry Tools
../include/mco/geometric/projective_geometry_utilities.h

# Parametric Integer Linear Programming
../include/mco/pilp/pilp_dual_benson.h
../include/mco/pilp/coin.h


#### Sources ####

# Benchmarks
benchmarks/lp_parser.cpp

# Benson Dual
generic/benson_dual/ove_cdd.cpp
generic/benson_dual/ove_node_lists.cpp
generic/benson_dual/ove_edge_lists.cpp
generic/benson_dual/ove_fp_v2.cpp
)

add_library(mco STATIC ${SOURCE_FILES})
