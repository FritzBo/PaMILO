include_directories(../include)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})

set(SOURCE_FILES

#### Header ####

# Core
../include/mco/basic/point.h
../include/mco/basic/pareto_point_comparator.h
../include/mco/basic/componentwise_point_comparator.h
../include/mco/basic/equality_point_comparator.h
../include/mco/basic/lex_point_comparator.h
../include/mco/basic/abstract_solver.h
../include/mco/basic/abstract_graph_instance.h
../include/mco/basic/weight_function_adaptors.h
../include/mco/basic/utility.h

# Assignment
../include/mco/ap/basic/abstract_ap_solver.h
../include/mco/ap/basic/ap_generator.h
../include/mco/ap/basic/ap_instance.h
../include/mco/ap/basic/lex_hungarian.h
../include/mco/ap/benson_dual/ap_benson_dual_solver.h
../include/mco/ap/benson_dual/ap_grb_benson_dual_solver.h
../include/mco/ap/brute_force/ap_brute_force_solver.h
../include/mco/ap/molp_solver/assignment_molp_model.h

# Benchmarks
../include/mco/benchmarks/abstract_objective_vector_generator.h
../include/mco/benchmarks/correlated_random.h
../include/mco/benchmarks/graph_generators.h
../include/mco/benchmarks/graph_parser.h
../include/mco/benchmarks/kurz_parser.h
../include/mco/benchmarks/mcap_parser.h
../include/mco/benchmarks/pge10_parser3_ap.h
../include/mco/benchmarks/random_objective_generator.h
../include/mco/benchmarks/temporary_graphs_parser.h

# MO Shortest Path
../include/mco/ep/basic/ep_instance.h
../include/mco/ep/basic/abstract_ep_solver.h
../include/mco/ep/brum_shier/ep_solver_bs.h
../include/mco/ep/brum_shier/ep_weighted_bs.h
../include/mco/ep/martins/label.h
../include/mco/ep/martins/martins_smart.h
../include/mco/ep/martins/martins.h
../include/mco/ep/martins/weighted_martins.h
../include/mco/ep/tsaggouris/ep_solver_tsaggouris_approx.h
../include/mco/ep/warburton/ep_solver_warburton_approx.h
../include/mco/ep/warburton/product.h
../include/mco/ep/basic/binary_heap.h
../include/mco/ep/basic/dijkstra.h
../include/mco/ep/dual_benson/ep_dual_benson.h


# MO Spanning Tree
../include/mco/est/basic/abstract_est_solver.h
../include/mco/est/basic/kruskal_st_solver.h
../include/mco/est/2tree/est_solver2_trees.h
../include/mco/est/dual_benson/est_dual_benson_scalarizer.h

# Benson Dual
../include/mco/generic/benson_dual/abstract_online_vertex_enumerator.h
../include/mco/generic/benson_dual/dual_benson_scalarizer.h
../include/mco/generic/benson_dual/ove_cdd.h
../include/mco/generic/benson_dual/ove_node_lists.h
../include/mco/generic/benson_dual/ove_edge_lists.h
../include/mco/generic/benson_dual/ove_fp_v2.h

# Geometry Tools
../include/mco/geometric/projective_geometry_utilities.h

# MO Linear Programming
../include/mco/molp/basic/molp_model.h
../include/mco/molp/benson/abstract_vertex_container.h
../include/mco/molp/benson/benson_exact.h
../include/mco/molp/benson/benson.h
../include/mco/molp/weightset/weightset_solver.h


#### Sources ####

# Assignment
ap/basic/ap_generator.cpp
ap/basic/lex_hungarian.cpp
ap/benson_dual/ap_benson_dual_solver.cpp
ap/benson_dual/ap_grb_benson_dual_solver.cpp
ap/brute_force/ap_brute_force_solver.cpp
ap/molp_solver/assignment_molp_model.cpp

# Benchmarks
benchmarks/correlated_random.cpp
benchmarks/graph_generators.cpp
benchmarks/graph_parser.cpp
benchmarks/kurz_parser.cpp
benchmarks/mcap_parser.cpp
benchmarks/pge10_parser3_ap.cpp
benchmarks/random_objective_generator.cpp
benchmarks/temporary_graphs_parser.cpp

# MO Shortest Path
ep/basic/ep_instance.cpp
ep/brum_shier/ep_solver_bs.cpp
ep/brum_shier/ep_weighted_bs.cpp
ep/martins/label.cpp
ep/martins/martins_smart.cpp
ep/martins/martins.cpp
ep/martins/weighted_martins.cpp
ep/tsaggouris/ep_solver_tsaggouris_approx.cpp
ep/warburton/ep_solver_warburton_approx.cpp
ep/basic/dijkstra.cpp

# MO Spanning Tree
est/basic/kruskal_st_solver.cpp
est/2tree/est_solver2_trees.cpp
est/dual_benson/est_dual_benson_scalarizer.cpp

# Benson Dual
generic/benson_dual/ove_cdd.cpp
generic/benson_dual/ove_node_lists.cpp
generic/benson_dual/ove_edge_lists.cpp
generic/benson_dual/ove_fp_v2.cpp

# MO Linear Programming
molp/basic/molp_model.cpp
molp/benson/benson_exact.cpp
molp/benson/benson.cpp
molp/weightset/weightset_solver.cpp
)

add_library(mco STATIC ${SOURCE_FILES})