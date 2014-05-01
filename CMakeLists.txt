cmake_minimum_required(VERSION 2.8)
PROJECT(SE3Tools)
cmake_policy(SET CMP0020 NEW)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" )
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
#find_package(Qt5 COMPONENTS Core Widgets REQUIRED)
#find_package(VTK REQUIRED)
#include(${VTK_USE_FILE})

add_subdirectory(SE3Tools)
add_subdirectory(SE3Tools_Test)

find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
			${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	add_custom_target(doc	${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
				WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
				COMMENT "Generating API Documentation with Doxygen" VERBATIM )
else(DOXYGEN_FOUND)
	MESSAGE(STATUS "Doxygen not found. Documentation will not be built.")
endif(DOXYGEN_FOUND)