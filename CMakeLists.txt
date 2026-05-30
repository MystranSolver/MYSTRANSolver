# CMakeLists.txt for compiling MYSTRAN with gfortran in a GNU environment
# based on an older CMakeLists.txt by ceanwang@gmail.com
# made to work again by Bruno Borges Paschoalinoto (2020)

# set up basic project info
cmake_minimum_required(VERSION 3.18)
include(CheckFunctionExists)

project(Mystran)
enable_language(Fortran)

if(WIN32)
  enable_language(RC)
endif()

# ! --- hijack mode for extern blas begin--- !
# ==============================================================================
# ABI SAFETY CHECK
# ==============================================================================
include(CheckFortranSourceCompiles)
check_fortran_source_compiles("
      program abi_check
      use iso_c_binding, only: c_int
      implicit none
      integer(c_int) :: n
      if (kind(n) /= 4) stop 1
      end program
" FORTRAN_C_INT_IS_4BYTE)

if(NOT FORTRAN_C_INT_IS_4BYTE)
  message(FATAL_ERROR
    \"\n*** BLAS ABI MISMATCH DETECTED ***\n\"
    \"Your Fortran compiler's C_INT is not 4 bytes. \"
    \"This means standard OpenBLAS (LP64) will corrupt the stack.\n\")
endif()
# ! --- hijack mode for extern blas end--- !

# ! --- hijack mode for extern blas begin--- !
# ==============================================================================
# FIND OPENBLAS
# ==============================================================================
# set(OpenBLAS_DIR "C:/gcc/openblas32/lib/cmake/openblas")
find_package(OpenBLAS QUIET)

if(OpenBLAS_FOUND)
    message(STATUS "Using OpenBLAS: ${OpenBLAS_LIBRARIES}")
    set(MYSTRAN_USE_OPENBLAS ON)
    set(enable_internal_blaslib OFF CACHE BOOL "" FORCE)
    set(TPL_BLAS_LIBRARIES "${OpenBLAS_LIBRARIES}" CACHE FILEPATH "" FORCE)
else()
    message(STATUS "OpenBLAS not found, will use internal BLAS")
    set(MYSTRAN_USE_OPENBLAS OFF)
    set(enable_internal_blaslib ON CACHE BOOL "" FORCE)
    unset(TPL_BLAS_LIBRARIES CACHE)
endif()
# ! --- hijack mode for extern blas end--- !

# basic compiler and output options
set(CMAKE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/Source")
set(BLAS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/BLAS")
set(PROJECT_BINARY_DIR "${PROJECT_SOURCE_DIR}/Binaries")
set(CMAKE_FLAGS "-c -fbacktrace")

# set(CMAKE_BUILD_TYPE Debug)

# set some dirs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/mod")
include_directories("${CMAKE_SOURCE_DIR}/INCLUDE")

# uncomment this to debug
# set(CMAKE_VERBOSE_MAKEFILE true)

# register the Profile build type, just in case
set(CMAKE_CONFIGURATION_TYPES
  Debug
  Release
  RelWithDebInfo
  MinSizeRel
  Profiling
  Deterministic
  CACHE STRING "" FORCE
)

# !--- ndebug release fix --- begin!
# Some Intel ifx Release setups conflict with forced NDEBUG definitions coming
# from default CMake release flags. Allow stripping NDEBUG for reproducible
# behavior across toolchains.
option(MYSTRAN_DISABLE_NDEBUG "Remove NDEBUG from CMake release-style flags" ON)
if(MYSTRAN_DISABLE_NDEBUG)
  foreach(lang C CXX Fortran)
    foreach(cfg RELEASE RELWITHDEBINFO MINSIZEREL)
      set(_var "CMAKE_${lang}_FLAGS_${cfg}")
      if(DEFINED ${_var})
        set(_flags "${${_var}}")
        string(REPLACE "-DNDEBUG" "" _flags "${_flags}")
        string(REPLACE "/DNDEBUG" "" _flags "${_flags}")
        string(REPLACE "-D NDEBUG" "" _flags "${_flags}")
        string(REPLACE "/D NDEBUG" "" _flags "${_flags}")
        string(REGEX REPLACE "[ \t]+" " " _flags "${_flags}")
        string(STRIP "${_flags}" _flags)
        set(${_var} "${_flags}" CACHE STRING "" FORCE)
      endif()
    endforeach()
  endforeach()
endif()
# !--- ndebug release fix --- end!

# suppress cmake warnings for superlu
if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
  set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")
endif()

# recommend the appropriate make command
if(WIN32)
  set(RECOMMENDED_MAKE "mingw32-make")
else()
  set(RECOMMENDED_MAKE "make")
endif()

# submodules (i.e. SuperLU)
find_package(Git QUIET)

if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  option(GIT_SUBMODULE "Check submodules during build" ON)

  if(GIT_SUBMODULE)
    message(STATUS "Updating submodules")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --force
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT
    )

    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(
        FATAL_ERROR
        "git submodule update --init --recursive --force failed"
      )
    endif()
  endif()
endif()

set(SUPERLU_DIR "${PROJECT_SOURCE_DIR}/superlu")
set(SUPERLU_MT_DIR "${PROJECT_SOURCE_DIR}/superlu_mt")
set(SUPERLU_PATCHES_DIR "${PROJECT_SOURCE_DIR}/superlu_mt_patches")

set(enable_examples OFF FORCE)
set(enable_tests OFF FORCE)
set(enable_doc OFF FORCE)

if(USE_SUPERLU_MT)
  message(STATUS "Will link against faster SuperLU_MT.")
  set(PLAT "_OPENMP" CACHE STRING "threading flavor _PTHREAD/_OPENMP" FORCE)
  message(STATUS "Applying patches to SuperLU_MT...")
  file(GLOB_RECURSE PATCH_FILES RELATIVE "${SUPERLU_PATCHES_DIR}" "${SUPERLU_PATCHES_DIR}/*")

  foreach(file IN LISTS PATCH_FILES)
    set(src "${SUPERLU_PATCHES_DIR}/${file}")
    set(dst "${SUPERLU_MT_DIR}/${file}")

    # ensure the destination directory exists
    get_filename_component(dst_dir "${dst}" DIRECTORY)
    file(MAKE_DIRECTORY "${dst_dir}")

    # print the file being copied
    message(STATUS "  Patching: ${file}")

    # copy the file
    file(COPY "${src}" DESTINATION "${dst_dir}")
  endforeach()

  message(STATUS "SuperLU_MT patched.")

  if(NOT EXISTS "${SUPERLU_MT_DIR}/CMakeLists.txt")
    message(
      FATAL_ERROR
      "The SuperLU_MT submodule was not downloaded. You have to do it manually!"
    )
  endif()

  if(WIN32)
    # force static linking to OpenMP library
    # temporarily set CMAKE_FIND_LIBRARY_SUFFIXES to prefer .a over .dll.a
    set(CMAKE_FIND_LIBRARY_SUFFIXES_BACKUP ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

    set(OpenMP_C_FLAGS "-fopenmp" CACHE STRING "" FORCE)
    set(OpenMP_C_LIB_NAMES "gomp" CACHE STRING "" FORCE)

    # find the static gomp library explicitly
    find_library(GOMP_STATIC_LIB NAMES libgomp.a gomp PATHS /mingw64/lib NO_DEFAULT_PATH)

    if(GOMP_STATIC_LIB)
      set(OpenMP_gomp_LIBRARY "${GOMP_STATIC_LIB}" CACHE FILEPATH "" FORCE)
    endif()

    find_package(OpenMP REQUIRED)

    # restore original library suffixes
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_BACKUP})
  endif()

  include_directories("${SUPERLU_MT_DIR}/SRC")
  add_subdirectory(${SUPERLU_MT_DIR})
  set(SLU_DRIVER "${SUPERLU_MT_DIR}/SRC/c_fortran_pdgssv.c")
else()
  message(STATUS "Will link against regular SuperLU.")

  if(NOT EXISTS "${SUPERLU_DIR}/CMakeLists.txt")
    message(
      FATAL_ERROR
      "The SuperLU submodule was not downloaded. You have to do it manually!"
    )
  endif()

  include_directories("${SUPERLU_DIR}/SRC")
  add_subdirectory(${SUPERLU_DIR})
  set(SLU_DRIVER "${SUPERLU_DIR}/FORTRAN/c_fortran_dgssv.c")
endif()

# ! --- hijack mode for extern blas begin--- !
if(NOT MYSTRAN_USE_OPENBLAS)
# f2c stuff
set(F2C_DIR "${PROJECT_SOURCE_DIR}/f2c")
set(F2C_INCLUDE_DIR "${F2C_DIR}/include")
set(F2C_FN "${F2C_DIR}/libf2c.zip")
set(F2C_URL "https://www.netlib.org/f2c/libf2c.zip")

# download f2c
if(NOT EXISTS ${F2C_DIR})
  message(STATUS "Downloading libf2c source from ${F2C_URL}...")
  make_directory("${F2C_DIR}")
  file(DOWNLOAD ${F2C_URL} ${F2C_FN} TIMEOUT 60 STATUS DOWNLOAD_STATUS)

  # Check if download was successful.
  list(GET DOWNLOAD_STATUS 0 STATUS_CODE)
  list(GET DOWNLOAD_STATUS 1 ERROR_MESSAGE)

  if(${STATUS_CODE} EQUAL 0)
    message(STATUS "Done downloading libf2c.")
  else()
    # Exit CMake if the download failed, printing the error message.
    file(REMOVE_RECURSE ${F2C_DIR})
    message(FATAL_ERROR "Error downloading libf2c: ${ERROR_MESSAGE}")
  endif()
endif()

# extract libf2c source
file(ARCHIVE_EXTRACT INPUT ${F2C_FN} DESTINATION ${F2C_DIR})

# prepare libf2c header files
file(GLOB_RECURSE F2C_PREHEADERS "${F2C_DIR}/*.h0")

foreach(H0 ${F2C_PREHEADERS})
  string(REGEX REPLACE "[.]h0$" ".h" H0_R ${H0})
  file(RENAME "${H0}" "${H0_R}")
  file(COPY "${H0_R}" DESTINATION "${F2C_INCLUDE_DIR}")
endforeach()

# get a load of this: f2c generates its own "arith.h" on the fly
# so we gotta compile arithchk and run it
set(F2C_ARITHCHK_SRC "${F2C_DIR}/arithchk.c")
set(F2C_ARITHCHK_BIN "${F2C_DIR}/arithchk")

if(WIN32)
  set(F2C_ARITHCHK_BIN "${F2C_ARITHCHK_BIN}.exe")
endif()

set(F2C_ARITH_H "${F2C_INCLUDE_DIR}/arith.h")
# --- feast_add --- begin !
set_source_files_properties(
  ${F2C_ARITHCHK_SRC} PROPERTIES COMPILE_FLAGS "-DNO_LONG_LONG -DNO_FPINIT"
)
add_executable(arithchk ${F2C_ARITHCHK_SRC})
target_link_libraries(arithchk m)
set_target_properties(
  arithchk PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${F2C_DIR}
)
add_custom_command(
  OUTPUT ${F2C_ARITH_H}
  COMMAND ${F2C_ARITHCHK_BIN} > ${F2C_ARITH_H}
  DEPENDS ${F2C_ARITHCHK_BIN}
)

# add libf2c to the compilation procedures
include_directories(${F2C_INCLUDE_DIR})
file(GLOB_RECURSE F2C_CFILES "${F2C_DIR}/*.c")
add_definitions(-DINTEGER_STAR_8)
add_library(f2c ${F2C_CFILES} ${F2C_ARITH_H})

# add some extra win32 flags for libf2c
if(WIN32)
  add_definitions(-DUSE_CLOCK -DMSDOS)
endif()
endif() # NOT MYSTRAN_USE_OPENBLAS
# ! --- hijack mode for extern blas end--- !

# set some extra vars for MSYS builds to make the binary portable
if(WIN32)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
endif()

# collect all fortran source files
file(GLOB_RECURSE ALL_FORTRAN_FILES
  "${CMAKE_SOURCE_DIR}/*.f"
  "${CMAKE_SOURCE_DIR}/*.F"
  "${CMAKE_SOURCE_DIR}/*.f90"
  "${CMAKE_SOURCE_DIR}/*.F90"
  "${CMAKE_SOURCE_DIR}/*.f95"
  "${CMAKE_SOURCE_DIR}/*.F95"
  "${CMAKE_SOURCE_DIR}/*.f03"
  "${CMAKE_SOURCE_DIR}/*.F03"
)

set_source_files_properties(
  "${PROJECT_SOURCE_DIR}/Source/LK4/EIGRL_EXTRACT_SOLVERS.F90"
  PROPERTIES Fortran_PREPROCESS ON COMPILE_OPTIONS "-cpp"
)


option(MYSTRAN_USE_EXTERNAL_FEAST "Enable native FEAST backend on condensed dense problem" OFF)
set(MYSTRAN_FEAST_EXTRA_LIBS "" CACHE STRING "Additional libraries for native FEAST backend")
# --- feast_add --- end !

# ! --- hijack mode for extern blas begin--- !
# BLAS configuration
if(MYSTRAN_USE_OPENBLAS)
    message(STATUS "Using OpenBLAS: ${TPL_BLAS_LIBRARIES}")
    set(missing_blas_src "${BLAS_SOURCE_DIR}/XERBLA.f")
else()
    # same BLAS-finding subroutine as SuperLU
    if(NOT enable_internal_blaslib)
      if(TPL_BLAS_LIBRARIES)
        set(BLAS_FOUND TRUE)
      else()
        find_package(BLAS)

        if(BLAS_FOUND)
          set(TPL_BLAS_LIBRARIES "${BLAS_LIBRARIES}" CACHE FILEPATH
            "Set from FindBLAS.cmake BLAS_LIBRARIES." FORCE)
        endif()
      endif()
    endif()

    if(NOT BLAS_FOUND)
      message(STATUS "BLAS not found, building local BLAS")
      list(
        APPEND blas_fns dgemm dgemv dlamch dlanst dscal dsteqr dsterf dswap dtrsm
        dtrtri ilaenv lsame xerbla
      )

      foreach(fname IN LISTS blas_fns)
        # message(STATUS "Checking for BLAS subr ${fname}")
        check_function_exists("${fname}" BLAS_FN_EXISTS)
        set(BLAS_FN_EXISTS CACHE "1" STRING)

        if(NOT BLAS_FN_EXISTS)
          # message(STATUS "BLAS subr ${fname} not found")
          string(TOUPPER ${fname} fname_upper)
          list(APPEND missing_blas_src "${BLAS_SOURCE_DIR}/${fname_upper}.f")
          list(APPEND missing_blas_fns ${fname})
        endif()

        unset(BLAS_FN_EXISTS CACHE)
      endforeach()

      # if any subroutines have bene found, create an inner blas library
      list(LENGTH missing_blas_fns MISSING_FNS_TOTAL)

      if(MISSING_FNS_TOTAL GREATER 0)
        if(MISSING_FNS_TOTAL GREATER 1)
          message(
            STATUS
            "BLAS subrs (${missing_blas_fns}) are absent and will be built locally."
          )
        else()
          message(
            STATUS
            "BLAS subr ${missing_blas_fns} is absent and will be built locally."
          )
        endif()
      endif()
    else()
      message(STATUS "Using system BLAS (${TPL_BLAS_LIBRARIES})")
    endif()
endif()
# ! --- hijack mode for extern blas end--- !

# prepare the main executable, linked against the specifics and the m
# it appears utils used to be a module, but that is no longer the case?
# file(GLOB UTIL_FILES "${CMAKE_SOURCE_DIR}/UTIL/*.f*")
# file(GLOB MAIN_FILES "${CMAKE_SOURCE_DIR}/MAIN/*.[fF]*")
add_executable(
  mystran
  ${ALL_FORTRAN_FILES}
  ${SLU_DRIVER}
  ${missing_blas_src}
)

# determine which superlu variant to link against
if(USE_SUPERLU_MT)
  if(PLAT STREQUAL "_PTHREAD")
    message(WARNING "We recommend using OpenMP, not pthread!")
    if(MYSTRAN_USE_OPENBLAS)
      target_link_libraries(mystran superlu_mt_PTHREAD ${TPL_BLAS_LIBRARIES})
    else()
      target_link_libraries(mystran superlu_mt_PTHREAD f2c)
    endif()
  elseif(PLAT STREQUAL "_OPENMP")
    if(MYSTRAN_USE_OPENBLAS)
      target_link_libraries(mystran superlu_mt_OPENMP ${TPL_BLAS_LIBRARIES})
    else()
      target_link_libraries(mystran superlu_mt_OPENMP f2c)
    endif()
  endif()

  target_compile_definitions(mystran PRIVATE USE_SUPERLU_MT)
else()
  if(MYSTRAN_USE_OPENBLAS)
    target_link_libraries(mystran superlu ${TPL_BLAS_LIBRARIES})
  else()
    target_link_libraries(mystran superlu f2c)
  endif()
endif()


if(MYSTRAN_USE_EXTERNAL_FEAST)
  target_compile_definitions(mystran PRIVATE MYSTRAN_HAVE_EXTERNAL_FEAST)
  foreach(_feast_lib IN LISTS MYSTRAN_FEAST_EXTRA_LIBS)
    get_filename_component(_feast_dir "${_feast_lib}" DIRECTORY)
    if(EXISTS "${_feast_dir}")
      target_include_directories(mystran PRIVATE "${_feast_dir}")
    endif()
  endforeach()
  if(MYSTRAN_FEAST_EXTRA_LIBS)
    separate_arguments(MYSTRAN_FEAST_EXTRA_LIBS)
    target_link_libraries(mystran ${MYSTRAN_FEAST_EXTRA_LIBS})
  endif()
endif()

set_target_properties(
  mystran PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

if(WIN32)
  target_sources(mystran PRIVATE mystran.rc)
endif()

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_Fortran_FLAGS "-Wall -Wextra \
    -Wno-unused-variable -Wno-unused-label -Wno-unused-parameter -Wno-tabs \
    -Wno-compare-reals -Wno-character-truncation -Wno-unused-dummy-argument \
    -Wmaybe-uninitialized -Wrealloc-lhs -fbounds-check -fcheck=all")
  set(CMAKE_Fortran_FLAGS_DEBUG
    "-g -O0 -fcheck=all -fbacktrace -fbounds-check \
    -fno-inline -fno-ipa-sra -fno-ipa-cp -fno-optimize-sibling-calls")
  set(CMAKE_Fortran_FLAGS_DETERMINISTIC
    "-O0 -g \
     -fno-fast-math \
     -ffp-contract=off \
     -fno-unsafe-math-optimizations \
     -fno-associative-math \
     -fno-reciprocal-math \
     -frounding-math"
  )
  set(CMAKE_C_FLAGS_DETERMINISTIC
    "-O0 -g \
     -fno-fast-math \
     -ffp-contract=off \
     -frounding-math"
  )
  set(CMAKE_CXX_FLAGS_DETERMINISTIC
    "${CMAKE_C_FLAGS_DETERMINISTIC}"
  )
  set(CMAKE_PROFILING_FLAGS "-O2 -g -fno-inline -fno-omit-frame-pointer -fno-inline-functions")
  set(CMAKE_C_FLAGS_PROFILING "${CMAKE_PROFILING_FLAGS}")
  set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_PROFILING_FLAGS}")
  set(CMAKE_Fortran_FLAGS_PROFILING "${CMAKE_PROFILING_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS_PROFILING "")
  set(CMAKE_SHARED_LINKER_FLAGS_PROFILING "")
endif()

# issue a couple messages about compilation
include(ProcessorCount)
ProcessorCount(NCPU)
message(STATUS "You can now compile MYSTRAN with `cmake --build .`")

if(NOT NCPU EQUAL 0)
  message(STATUS "Compile faster by passing the -j${NCPU} flag.")
endif()
