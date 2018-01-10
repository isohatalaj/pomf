# Find the uthash header

include(LibFindMacros)

libfind_pkg_check_modules(uthash_PKGCONF)

find_path(uthash_INCLUDE_DIR
  NAMES uthash.h
  PATHS ${uthash_PKGCONF_INCLUDE_DIRS}
)

set(cxsparse_PROCESS_INCLUDES uthash_INCLUDE_DIR)
libfind_process(uthash)

