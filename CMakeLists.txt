# CMake 最低版本要求
cmake_minimum_required(VERSION 3.11...3.16)

# 项目名称和设置
project(FEM_IPC
    VERSION 0.1
    DESCRIPTION "Hydrostatic Object Project"
    LANGUAGES CXX)

# 选项设置
OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()    

# 设置Windows的导出符号
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# 定义外部项目目录
set(EXTERN_PROJECTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern)


# 设置包含路径
list (APPEND EXTRA_INCLUDES 
    ${EXTERN_PROJECTS_DIR}/eigen
    ${EXTERN_PROJECTS_DIR}/cppAD/include  # 指向cppAD的include目录
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Energy
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FormatMesh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/IPC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Objects
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Simulator
    ${CMAKE_CURRENT_SOURCE_DIR}/include/MLS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/Utils
)

# 添加可执行文件及其源文件和头文件
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
    src/Energy/ExternalEnergy.cpp
    src/Energy/BarrierEnergy.cpp
    src/IPC/distance.cpp
    src/IPC/CCD.cpp
    src/FormatMesh/objMesh.cpp
    src/MLS/MLS.cpp
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
)

# 组织源文件和头文件到Visual Studio的项目组
source_group("Source Files\\Utils" FILES src/Utils/utils.cpp src/Utils/tools.cpp)
source_group("Source Files\\Objects" FILES src/Objects/mesh.cpp)
source_group("Source Files\\Energy" FILES 
    src/Energy/InertiaEnergy.cpp 
    src/Energy/ElasticEnergy.cpp
    src/Energy/ExternalEnergy.cpp
    src/Energy/BarrierEnergy.cpp)
source_group("Source Files\\Simulator" FILES 
    src/Simulator/simulator.cpp 
    src/Simulator/simulator_ABD.cpp)
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

# 针对Windows的特定编译选项
if(WIN32)
    target_compile_options(${PROJECT_NAME} PRIVATE /bigobj)
endif()

# 使用C++14标准
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_14)

# 添加包含路径
target_include_directories(${PROJECT_NAME} PRIVATE ${EXTRA_INCLUDES})

# 链接库，如果有的话
target_link_libraries (${PROJECT_NAME} ${EXTRA_LIBS}) 

# 定义预处理宏
target_compile_definitions(${PROJECT_NAME} PUBLIC -DROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")