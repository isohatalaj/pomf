cmake_minimum_required(VERSION 2.8.4)
project(pomf)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# find_package(LAPACK) may fail on MinGW. In this case, give
# the path to blas and lapack when invoking cmake, i.e.
#
#   cmake -DLAPACK_LIBRARIES=<full-path-to/liblapack.a> -DBLAS_LIBRARIES=...
#
# Further, the fortran runtime library that was used to build blas and lapack 
# may need to be given. If linker complains e.g. about undefined _gfortran_...
# references, add -DFORTRAN_RT_LIBRARY=<full-path-to/your-fortran-runtime.a>.
# On MinGW, one may thus need, for example:
#
#   cmake -G "MinGW Makefiles" -DLAPACK_LIBRARIES=c:/MinGW/lib/liblapack.a  \
#      -DBLAS_LIBRARIES=c:/MinGW/lib/libblas.a \
#      -DFORTRAN_RT_LIBRARY=c:/MinGW/lib/gcc/mingw32/6.3.0/libgfortran.a;c:/MinGW/lib/gcc/mingw32/6.3.0/libquadmath.a
#

find_package(GSL REQUIRED)
find_package(uthash REQUIRED)
find_package(LAPACK REQUIRED)
find_package(cxsparse)
find_package(umfpack)

if (cxsparse_FOUND AND umfpack_FOUND)
  set(sparselibs_FOUND TRUE)
  set(sparselibs_INCLUDE_DIRS ${cxsparse_INCLUDE_DIRS} ${umfpack_INCLUDE_DIRS})
  set(sparselibs_LIBRARIES ${cxsparse_LIBRARIES} ${umfpack_LIBRARIES})
else (cxsparse_FOUND AND umfpack_FOUND)
  set(sparselibs_FOUND FALSE)
endif (cxsparse_FOUND AND umfpack_FOUND)

if (DEFINED FORTRAN_RT_LIBRARY)
  set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${FORTRAN_RT_LIBRARY})
  set(BLAS_LIBRARIES ${BLAS_LIBRARIES} ${FORTRAN_RT_LIBRARY})	
endif (DEFINED FORTRAN_RT_LIBRARY)

include_directories(${GSL_INCLUDE_DIRS} ${uthash_INCLUDE_DIRS} ${sparselibs_INCLUDE_DIRS})

if (UNIX)
  find_library(MATH_LIBRARY m)
endif (UNIX)

if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O2 -Wall -Werror")
endif ()

add_subdirectory(src)
add_subdirectory(util)

