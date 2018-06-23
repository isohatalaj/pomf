
set(LIBSMUTIL_SOURCE util.h util.c rgrid.h rgrid.c spgrid.h spgrid.c)
set(LIBSMUTIL_INCLUDE util.h rgrid.h spgrid.h)

if (sparselibs_FOUND)
   set(LIBSMUTIL_SOURCE ${LIBSMUTIL_SOURCE} smaux.h smaux.c)
   set(LIBSMUTIL_INCLUDE ${LIBSMUTIL_INCLUDE} smaux.h)
endif (sparselibs_FOUND)

set(all_libs ${MATH_LIBRARY} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${sparselibs_LIBRARIES})

add_library(smutil SHARED ${LIBSMUTIL_SOURCE})
target_link_libraries(smutil ${all_libs})




