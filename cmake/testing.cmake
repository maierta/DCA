################################################################################
# Description.
#
# References: - https://github.com/ALPSCore/ALPSCore
################################################################################

include(CMakeParseArguments)

function(add_gtest)
  set(options ADVANCED MPI CUDA GTEST_MAIN)
  set(oneValueArgs NAME MPI_NUMPROC)
  set(multiValueArgs INCLUDES SOURCES LIBS)
  cmake_parse_arguments(ADD_GTEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if (ADD_GTEST_ADVANCED AND NOT DCA_TESTS_INCLUDE_ADVANCED)
    # Do not build this advanced test.
    return()
  endif()

  if (ADD_GTEST_MPI AND NOT DCA_HAVE_MPI)
    message(WARNING "Requested MPI for test ${ADD_GTEST_NAME} but MPI is not available. Test is skipped.")
    return()
  endif()

  if (ADD_GTEST_CUDA AND NOT DCA_HAVE_CUDA)
    message(WARNING "Requested CUDA for test ${ADD_GTEST_NAME} but CUDA is not available. Test is skipped.")
    return()
  endif()

  add_executable(${ADD_GTEST_NAME} ${ADD_GTEST_NAME}.cpp "${ADD_GTEST_SOURCES}")
  
  if (ADD_GTEST_GTEST_MAIN)
    # Use gtest main.
    target_link_libraries(${ADD_GTEST_NAME} gtest_main "${ADD_GTEST_LIBS}")
  else()
    # Test has its own main.
    target_link_libraries(${ADD_GTEST_NAME} gtest "${ADD_GTEST_LIBS}")
  endif()

  # TODO: Pass GPU routines and GPU libraries as argument to add_gtest function?
  if (ADD_GTEST_CUDA)
       target_link_libraries(${ADD_GTEST_NAME}
                             DCA_GPU_routines
                             ${DCA_GPU_LIBRARIES})
       cuda_add_cublas_to_target(${ADD_GTEST_NAME})
  endif()

  target_include_directories(${ADD_GTEST_NAME} PRIVATE "${ADD_GTEST_INCLUDES}")

  if (ADD_GTEST_MPI)
    if (NOT DEFINED ADD_GTEST_MPI_NUMPROC)
      # message("MPI_NUMPROC not specified. Setting it to 4.")
      set(ADD_GTEST_MPI_NUMPROC 4)
    endif()
    
    add_test(
      NAME    ${ADD_GTEST_NAME}
      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${ADD_GTEST_MPI_NUMPROC} ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:${ADD_GTEST_NAME}>" ${MPIEXEC_POSTFLAGS})

  else()
    add_test(
      NAME    ${ADD_GTEST_NAME}
      COMMAND "$<TARGET_FILE:${ADD_GTEST_NAME}>")

  endif()
endfunction()

option(DCA_TESTS                  "Build DCA++'s tests."                        OFF)
option(DCA_TESTS_INCLUDE_ADVANCED "Include time- and resource-consuming tests." OFF)

if (DCA_TESTS)
  enable_testing()

  add_subdirectory(${gtest_DIR} ${PROJECT_BINARY_DIR}/gtest)
  
  target_compile_options(gtest      PRIVATE "-w")
  target_compile_options(gtest_main PRIVATE "-w")
  
  include_directories(${gtest_SOURCE_DIR}/include)

  add_subdirectory(${PROJECT_SOURCE_DIR}/test/unit)
endif()
