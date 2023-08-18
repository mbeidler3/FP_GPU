cmake_minimum_required (VERSION 2.8.8)

project (fpgpu_proj Fortran)

#  Set directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/bin)

#  Fortran specific settings.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/modules)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/lib)

add_definitions (-DDOUBLE_PRECISION)

option (USE_ACC "Use OpenACC" OFF)
if (${USE_ACC})
   add_definitions (-DACC)
endif()

add_subdirectory (src)
