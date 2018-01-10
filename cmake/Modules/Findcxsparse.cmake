# Find the cxsparse library

include(LibFindMacros)

libfind_pkg_check_modules(cxsparse_PKGCONF)

find_path(cxsparse_INCLUDE_DIR
  NAMES cs.h
  PATHS ${cxsparse_PKGCONF_INCLUDE_DIRS}
)

find_library(cxsparse_LIBRARY
  NAMES cxsparse
  PATHS ${cxsparse_PKGCONF_LIBRARY_DIRS}
)

set(cxsparse_PROCESS_INCLUDES cxsparse_INCLUDE_DIR)
set(cxsparse_PROCESS_LIBS cxsparse_LIBRARY)
libfind_process(cxsparse)

