
include_directories(${PROJECT_SOURCE_DIR}/util)
link_directories(${PROJECT_BINARY_DIR}/util)

set(LIBPOMF_SOURCE pomf.c model.c mcsimu.c)
set(LIBPOMF_INCLUDE pomf.h pomf_model.h pomf_mcsimu.h)

if (sparselibs_FOUND)
   set(LIBPOMF_SOURCE ${LIBPOMF_SOURCE} smsimu.c umfsolve.c fps.c)
   set(LIBPOMF_INCLUDE ${LIBPOMF_SOURCE} pomf_smsimu.h fps.h)
endif (sparselibs_FOUND)

set(all_libs ${MATH_LIBRARY} ${GSL_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${sparselibs_LIBRARIES})

add_library(pomf SHARED ${LIBPOMF_SOURCE})
target_link_libraries(pomf smutil ${all_libs})

add_executable(mcrun mcrun.c)
target_link_libraries(mcrun pomf smutil ${all_libs})

add_executable(mcstep mcstep.c)
target_link_libraries(mcstep pomf smutil ${all_libs})

if (sparselibs_FOUND)
   add_executable(smrun smrun.c)
   target_link_libraries(smrun pomf smutil ${all_libs})

   add_executable(steptest steptest.c)
   target_link_libraries(steptest pomf smutil ${all_libs})

   add_executable(fpstest fpstest.c)
   target_link_libraries(fpstest pomf smutil ${all_libs})
endif (sparselibs_FOUND)

