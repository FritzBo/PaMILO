include_directories(../include)
include_directories(../tclap)
include_directories(${GUROBI_INCLUDE_PATH})
include_directories(${COIN_INCLUDE_PATH})
include_directories(${OGDF_INCLUDE_PATH})

set(SOURCE_FILES
main.cpp

#### Core
basic/modules.h
basic/modules.cpp

#### Modules
modules/ep_benson_module.h
modules/ep_benson_module.cpp
modules/ep_martins_module.h
modules/ep_martins_module.cpp
)

include_directories(${BOOST_PO_INCLUDE})

add_executable(mco_cli ${SOURCE_FILES})

target_link_libraries(mco_cli mco)
target_link_libraries(mco_cli debug ${OGDF-DBG} optimized ${OGDF})
target_link_libraries(mco_cli debug ${COIN-DBG} optimized ${COIN})
target_link_libraries(mco_cli ${CDD})
target_link_libraries(mco_cli pthread)