cmake_minimum_required(VERSION 2.6)
PROJECT(SE3Tools)
cmake_policy(SET CMP0020 NEW)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" )
find_package(Eigen3 REQUIRED)
#find_package(Qt5 COMPONENTS Core Widgets REQUIRED)
#find_package(VTK REQUIRED)
#include(${VTK_USE_FILE})

add_subdirectory(SE3Tools)
add_subdirectory(SE3Tools_Test)