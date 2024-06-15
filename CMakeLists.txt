# Works with 3.11 and tested through 3.15
cmake_minimum_required(VERSION 3.11...3.16)

# Project name and a few useful settings. Other commands can pick up the results
project(FEM_IPC
    VERSION 0.1
    DESCRIPTION "Hrdrostatic Object Project"
    LANGUAGES CXX)
   
    
OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()    
    
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

set(EXTERN_PROJECTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern)

add_subdirectory(${EXTERN_PROJECTS_DIR})

if(WIN32)
    list (APPEND EXTRA_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME} ${EXTERN_PROJECTS_DIR}/eigen)
else()
    list (APPEND EXTRA_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME} ${EXTERN_PROJECTS_DIR}/eigen)
endif()






add_executable(
    ${PROJECT_NAME} 
    src/main.cpp
	src/utils.cpp
	src/mesh.cpp
	src/InertiaEnergy.cpp
	src/ElasticEnergy.cpp
	src/simulator.cpp
	src/GravityEnergy.cpp
	include/FEM_IPC/GravityEnergy.h
	include/FEM_IPC/simulator.h
	include/FEM_IPC/ElasticEnergy.h
	include/FEM_IPC/InertiaEnergy.h
	include/FEM_IPC/utils.h
	include/FEM_IPC/materials.h
	include/FEM_IPC/mesh.h)
    
    

if(WIN32)
    target_compile_options(${PROJECT_NAME} PRIVATE /bigobj)
endif()
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_14)
target_include_directories(${PROJECT_NAME} PRIVATE ${EXTRA_INCLUDES})
target_link_libraries (${PROJECT_NAME} ${EXTRA_LIBS}) 
target_compile_definitions(${PROJECT_NAME} PUBLIC -DROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
