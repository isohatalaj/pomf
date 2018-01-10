cmake_minimum_required(VERSION 2.8.4)
project(pomf)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(GSL REQUIRED)
find_package(cxsparse REQUIRED)
find_package(uthash REQUIRED)

if (UNIX)
  find_library(MATH_LIBRARY m)
endif (UNIX)

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O2 -ffast-math -Wall -Werror")
endif()

add_subdirectory(src)
add_subdirectory(util)

