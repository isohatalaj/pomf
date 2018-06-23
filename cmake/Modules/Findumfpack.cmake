# Find the umfpack library

include(LibFindMacros)

libfind_pkg_check_modules(umfpack_PKGCONF)

find_path(umfpack_INCLUDE_DIR
  NAMES umfpack.h
  PATHS ${umfpack_PKGCONF_INCLUDE_DIRS}
)

find_library(umfpack_LIBRARY
  NAMES umfpack
  PATHS ${umfpack_PKGCONF_LIBRARY_DIRS}
)

set(umfpack_PROCESS_INCLUDES umfpack_INCLUDE_DIR)
set(umfpack_PROCESS_LIBS umfpack_LIBRARY)
libfind_process(umfpack)

