
set(LIBSMUTIL_SOURCE util.h util.c rgrid.h rgrid.c spgrid.h spgrid.c smaux.h smaux.c)
set(LIBSMUTIL_INCLUDE util.h rgrid.h spgrid.h smaux.h)

add_library(smutil ${LIBSMUTIL_SOURCE})

set(all_libs ${MATH_LIBRARY})


