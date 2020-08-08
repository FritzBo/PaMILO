include_directories(../include)
include_directories(../tclap)
include_directories(${CPLEX_INCLUDE_DIRS})
include_directories(${CONCERT_INCLUDE_DIRS})
include_directories(${PROJECT_BINARY_DIR})

if(USE_CDD)
	include_directories(${CDD_INCLUDE_PATH})
endif()

set(SOURCE_FILES
main.cpp

#### Core
basic/modules.h
basic/modules.cpp

#### Modules
modules/pilp_benson_module.h
modules/pilp_benson_module.cpp
)

add_executable(pamilo_cli ${SOURCE_FILES})

target_link_libraries(pamilo_cli pamilo)
target_link_libraries(pamilo_cli ${ILOCPLEX_LIB})
target_link_libraries(pamilo_cli ${CONCERT_LIB})
target_link_libraries(pamilo_cli ${CPLEX_LIB})
if(USE_CDD)
	target_link_libraries(pamilo_cli ${CDD})
endif()
target_link_libraries(pamilo_cli pthread)
target_link_libraries(pamilo_cli dl)
