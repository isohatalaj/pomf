
include_directories(${PROJECT_SOURCE_DIR}/util)
link_directories(${PROJECT_BINARY_DIR}/util)

set(LIBPOMF_SOURCE pomf.c model.c mcsimu.c smsimu.c)
set(LIBPOMF_INCLUDE pomf.h pomf_model.h pomf_mcsimu.h pomf_smsimu.h)

add_library(pomf ${LIBPOMF_SOURCE})

set(all_libs ${MATH_LIBRARY} ${GSL_LIBRARIES} ${cxsparse_LIBRARIES})

add_executable(mcrun mcrun.c)
target_link_libraries(mcrun pomf smutil ${all_libs})

add_executable(smrun smrun.c)
target_link_libraries(smrun pomf smutil ${all_libs})


