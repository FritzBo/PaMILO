include_directories(../include)
include_directories(../tclap)
include_directories(${CPLEX_INCLUDE_DIRS})
include_directories(${CONCERT_INCLUDE_DIRS})

set(SOURCE_FILES
main.cpp

#### Core
basic/modules.h
basic/modules.cpp

#### Modules
modules/pilp_benson_module.h
modules/pilp_benson_module.cpp
)

#include_directories(${BOOST_PO_INCLUDE})

add_executable(mco_cli ${SOURCE_FILES})

target_link_libraries(mco_cli mco)
target_link_libraries(mco_cli ${ILOCPLEX_LIB_PATH})
target_link_libraries(mco_cli ${CONCERT_LIB_PATH})
target_link_libraries(mco_cli ${CPLEX_LIB_PATH})
target_link_libraries(mco_cli ${CDD})
#target_link_libraries(mco_cli m)
target_link_libraries(mco_cli pthread)
target_link_libraries(mco_cli dl)
