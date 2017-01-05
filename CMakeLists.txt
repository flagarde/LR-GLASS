########################################################
# cmake file for building LRGRPC
# @author Lagarde François, IPNL
########################################################
cmake_minimum_required(VERSION 2.8)
if(CMAKE_VERSION VERSION_LESS "2.8.12")
   add_definitions("-std=c++14")
else()
   add_compile_options(-std=c++14)
endif()

project( LRGRPC )
set( ${PROJECT_NAME}_VERSION_MAJOR 7 )
set( ${PROJECT_NAME}_VERSION_MINOR 0 )
set( ${PROJECT_NAME}_VERSION_PATCH 0 )

add_definitions( "-D_GLIBCXX_USE_CXX11_ABI=0 -ggdb -Wall -ansi -fno-strict-aliasing -fPIC -O3" )
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/source/ROOT/)
list( APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake/)
configure_file ("./source/include/Version.h.in" "${CMAKE_SOURCE_DIR}/build/include/Version.h")
include_directories(${CMAKE_SOURCE_DIR}/source/src/)
include_directories(${CMAKE_SOURCE_DIR}/source/include/)
include_directories(${CMAKE_SOURCE_DIR}/build/include/)

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})

AUX_SOURCE_DIRECTORY( ./source/src/ library_sources )
LIST( REMOVE_ITEM library_sources ./source/src/main.cpp)
ADD_LIBRARY( LIBLRGRPC STATIC ${library_sources} )
ADD_EXECUTABLE( LRGRPC ./source/src/main.cpp)
TARGET_LINK_LIBRARIES( LRGRPC LIBLRGRPC ${ROOT_LIBRARIES})
INSTALL( TARGETS LRGRPC DESTINATION ${CMAKE_SOURCE_DIR}/ )
