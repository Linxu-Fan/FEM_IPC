# CMake 最低版本要求
cmake_minimum_required(VERSION 3.11...3.16)

# 项目名称和设置
project(FEM_IPC
    VERSION 0.1
    DESCRIPTION "Hydrostatic Object Project"
    LANGUAGES CXX)


if(MSVC)
    add_compile_options("/utf-8")
    # 或者
    add_definitions(-D_UNICODE -DUNICODE)
endif()



# Enable C++17 (or the version CGAL requires)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# 选项设置
OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()



# Find CGAL package
find_package(CGAL CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)
find_package(OpenVDB CONFIG REQUIRED)


# 设置Windows的导出符号
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# 定义外部项目目录
set(EXTERN_PROJECTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern)

# 设置包含路径
list (APPEND EXTRA_INCLUDES
	${EXTERN_PROJECTS_DIR}/libigl/include 
    ${EXTERN_PROJECTS_DIR}/eigen
	${EXTERN_PROJECTS_DIR}/voro++/src
    ${EXTERN_PROJECTS_DIR}/cppAD/include  # 指向cppAD的include目录
	${VORO_PLUS_PLUS_INCLUDE_DIRS} 
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Energy
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FormatMesh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/IPC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Objects
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Simulator
    ${CMAKE_CURRENT_SOURCE_DIR}/include/MLS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Utils
	${CMAKE_CURRENT_SOURCE_DIR}/include/Crack
)


message(STATUS "TBB_INCLUDE_DIRS: ${TBB_INCLUDE_DIRS}")
message(STATUS "TBB_LIBRARIES: ${TBB_LIBRARIES}")



# === 添加可执行文件及其源文件和头文件 ===
add_executable(
    ${PROJECT_NAME} 
    src/main.cpp
    src/Utils/utils.cpp
    src/Utils/tools.cpp
    src/Objects/mesh.cpp
    src/Energy/InertiaEnergy.cpp
    src/Energy/ElasticEnergy.cpp
    src/Simulator/simulator.cpp
    src/Simulator/simulator_ABD.cpp
	src/Simulator/simulator_ABD_triMesh.cpp
    src/Energy/ExternalEnergy.cpp
    src/Energy/BarrierEnergy.cpp
    src/IPC/distance.cpp
    src/IPC/CCD.cpp
    src/FormatMesh/objMesh.cpp
    src/MLS/MLS.cpp
	src/Crack/extractCrack.cpp
    include/MLS/MLS.h
    include/Utils/tools.h
    include/Utils/utils.h
    include/FormatMesh/objMesh.h
    include/Simulator/simulator.h
    include/IPC/distance.h
    include/IPC/CCD.h
    include/Energy/BarrierEnergy.h
    include/Energy/ExternalEnergy.h
    include/Energy/ElasticEnergy.h
    include/Energy/InertiaEnergy.h	
    include/Objects/materials.h
    include/Objects/mesh.h
	include/Crack/extractCrack.h
)

# 添加 libigl 子目录
set(LIBIGL_WITH_CGAL ON CACHE BOOL "Use CGAL for mesh operations")
set(LIBIGL_WITH_OPENGL ON CACHE BOOL "Use OpenGL")
set(LIBIGL_WITH_OPENGL_GLFW ON CACHE BOOL "Use OpenGL with GLFW")
set(LIBIGL_WITH_OPENGL_GLFW_IMGUI OFF CACHE BOOL "Use ImGui")



add_subdirectory(${EXTERN_PROJECTS_DIR}/libigl)

# 组织源文件和头文件到 Visual Studio 的项目组
source_group("Source Files\\Utils" FILES src/Utils/utils.cpp src/Utils/tools.cpp)
source_group("Source Files\\Objects" FILES src/Objects/mesh.cpp)
source_group("Source Files\\Crack" FILES src/Crack/extractCrack.cpp)
source_group("Source Files\\Energy" FILES 
    src/Energy/InertiaEnergy.cpp 
    src/Energy/ElasticEnergy.cpp
    src/Energy/ExternalEnergy.cpp
    src/Energy/BarrierEnergy.cpp)
source_group("Source Files\\Simulator" FILES 
    src/Simulator/simulator.cpp 
    src/Simulator/simulator_ABD.cpp
	src/Simulator/simulator_ABD_triMesh.cpp)
source_group("Source Files\\IPC" FILES src/IPC/distance.cpp src/IPC/CCD.cpp)
source_group("Source Files\\FormatMesh" FILES src/FormatMesh/objMesh.cpp)
source_group("Source Files\\MLS" FILES src/MLS/MLS.cpp)

source_group("Header Files\\Utils" FILES include/Utils/tools.h include/Utils/utils.h)
source_group("Header Files\\FormatMesh" FILES include/FormatMesh/objMesh.h)
source_group("Header Files\\Simulator" FILES include/Simulator/simulator.h)
source_group("Header Files\\IPC" FILES include/IPC/distance.h include/IPC/CCD.h)
source_group("Header Files\\Energy" FILES 
    include/Energy/BarrierEnergy.h 
    include/Energy/ExternalEnergy.h 
    include/Energy/ElasticEnergy.h 
    include/Energy/InertiaEnergy.h)
source_group("Header Files\\Objects" FILES include/Objects/materials.h include/Objects/mesh.h)
source_group("Header Files\\MLS" FILES include/MLS/MLS.h)
source_group("Header Files\\Crack" FILES include/Crack/extractCrack.h)

# 针对 Windows 的特定编译选项
if(WIN32)
    target_compile_options(${PROJECT_NAME} PRIVATE /bigobj)
endif()

# 使用 C++14 标准
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_14)

# 添加包含路径
target_include_directories(${PROJECT_NAME} PRIVATE ${EXTRA_INCLUDES})


target_link_libraries(${PROJECT_NAME} PRIVATE TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenVDB::openvdb)


# 链接库，包括 libigl 和 TBB
target_link_libraries(${PROJECT_NAME} PUBLIC 
    ${EXTRA_LIBS}
    igl::core
	CGAL::CGAL
)

# 定义预处理宏
target_compile_definitions(${PROJECT_NAME} PUBLIC -DROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")