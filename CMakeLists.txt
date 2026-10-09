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

string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" _MYSTRAN_CPU_ARCH)
set(_MYSTRAN_DEFAULT_MARCH "")
set(_MYSTRAN_DEFAULT_MTUNE "")

if(_MYSTRAN_CPU_ARCH MATCHES "^(x86_64|amd64|x64)$")
  set(_MYSTRAN_DEFAULT_MARCH "x86-64")
  set(_MYSTRAN_DEFAULT_MTUNE "generic")
endif()

set(MYSTRAN_MARCH "${_MYSTRAN_DEFAULT_MARCH}" CACHE STRING "CPU instruction set (-march); empty for compiler default")
set(MYSTRAN_MTUNE "${_MYSTRAN_DEFAULT_MTUNE}" CACHE STRING "CPU tuning (-mtune); empty for compiler default")

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
include(CheckFortranCompilerFlag)
set(_MYSTRAN_CPU_FLAGS "")

foreach(_cpu_option MARCH MTUNE)
  if(NOT "${MYSTRAN_${_cpu_option}}" STREQUAL "")
    if(NOT "${MYSTRAN_${_cpu_option}}" MATCHES "^[A-Za-z0-9_.+-]+$")
      message(FATAL_ERROR "MYSTRAN_${_cpu_option} must be a single CPU name or empty.")
    endif()

    string(TOLOWER "${_cpu_option}" _cpu_option_lower)
    set(_cpu_flag "-${_cpu_option_lower}=${MYSTRAN_${_cpu_option}}")
    string(MD5 _cpu_flag_hash "${_cpu_flag}")

    foreach(_cpu_language C CXX Fortran)
      set(_cpu_check "_MYSTRAN_${_cpu_language}_CPU_${_cpu_flag_hash}")

      if(_cpu_language STREQUAL "C")
        check_c_compiler_flag("${_cpu_flag}" ${_cpu_check})
      elseif(_cpu_language STREQUAL "CXX")
        check_cxx_compiler_flag("${_cpu_flag}" ${_cpu_check})
      else()
        check_fortran_compiler_flag("${_cpu_flag}" ${_cpu_check})
      endif()

      if(NOT ${_cpu_check})
        message(FATAL_ERROR
          "${CMAKE_${_cpu_language}_COMPILER_ID} ${_cpu_language} compiler does not support ${_cpu_flag}. "
          "Choose a supported MYSTRAN_${_cpu_option}, or set it to an empty string for the compiler default.")
      endif()
    endforeach()

    list(APPEND _MYSTRAN_CPU_FLAGS "${_cpu_flag}")
    add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX,Fortran>:${_cpu_flag}>")
  endif()
endforeach()

message(STATUS "MYSTRAN CPU flags: ${_MYSTRAN_CPU_FLAGS}")

# basic compiler and output options
set(CMAKE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/Source")
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

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# suppress cmake warnings for superlu
if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
  set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")
endif()

# build static libraries throughout this tree
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libraries" FORCE)

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

# All git submodules and their patch overlays live under submodules/.
set(SUBMODULES_DIR "${PROJECT_SOURCE_DIR}/submodules")
set(SUPERLU_DIR "${SUBMODULES_DIR}/superlu")
set(SUPERLU_PATCHES_DIR "${SUBMODULES_DIR}/superlu_patches")
set(SUPERLU_MT_DIR "${SUBMODULES_DIR}/superlu_mt")
set(SUPERLU_MT_PATCHES_DIR "${SUBMODULES_DIR}/superlu_mt_patches")
set(METIS_DIR "${SUBMODULES_DIR}/metis")
set(GKLIB_DIR "${SUBMODULES_DIR}/GKlib")
set(METIS_PATCHES_DIR "${SUBMODULES_DIR}/metis_patches")
set(GKLIB_PATCHES_DIR "${SUBMODULES_DIR}/gklib_patches")
set(MYSTRAN_LAPACK_DIR "${SUBMODULES_DIR}/lapack")

if(NOT DEFINED TPL_ENABLE_METISLIB)
  set(TPL_ENABLE_METISLIB ON CACHE BOOL "Enable METIS for SuperLU")
endif()

if(TPL_ENABLE_METISLIB)
  if(NOT EXISTS "${GKLIB_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "The GKlib submodule was not downloaded.")
  endif()

  if(NOT EXISTS "${METIS_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "The METIS submodule was not downloaded.")
  endif()

  message(STATUS "Applying patches to GKlib...")
  file(GLOB_RECURSE GKLIB_PATCH_FILES RELATIVE "${GKLIB_PATCHES_DIR}" "${GKLIB_PATCHES_DIR}/*")

  foreach(file IN LISTS GKLIB_PATCH_FILES)
    set(src "${GKLIB_PATCHES_DIR}/${file}")
    set(dst "${GKLIB_DIR}/${file}")

    get_filename_component(dst_dir "${dst}" DIRECTORY)
    file(MAKE_DIRECTORY "${dst_dir}")

    message(STATUS "  Patching: ${file}")
    file(COPY "${src}" DESTINATION "${dst_dir}")
  endforeach()

  message(STATUS "GKlib patched.")

  message(STATUS "Applying patches to METIS...")
  file(GLOB_RECURSE METIS_PATCH_FILES RELATIVE "${METIS_PATCHES_DIR}" "${METIS_PATCHES_DIR}/*")

  foreach(file IN LISTS METIS_PATCH_FILES)
    set(src "${METIS_PATCHES_DIR}/${file}")
    set(dst "${METIS_DIR}/${file}")

    get_filename_component(dst_dir "${dst}" DIRECTORY)
    file(MAKE_DIRECTORY "${dst_dir}")

    message(STATUS "  Patching: ${file}")
    file(COPY "${src}" DESTINATION "${dst_dir}")
  endforeach()

  message(STATUS "METIS patched.")

  set(SHARED OFF CACHE BOOL "Build vendored GKlib/METIS as static libraries" FORCE)
  set(GKLIB_BUILD_APPS OFF CACHE BOOL "Build GKlib applications" FORCE)
  set(METIS_BUILD_PROGRAMS OFF CACHE BOOL "Build METIS command line programs" FORCE)
  set(GKLIB_PATH "${GKLIB_DIR}" CACHE PATH "GKlib source directory" FORCE)

  # Use a binary directory that is disjoint from the source tree. The MYSTRAN
  # build is performed in-source (CMAKE_BINARY_DIR == CMAKE_SOURCE_DIR), so
  # using "${CMAKE_BINARY_DIR}/metis/METIS" would make METIS' binary dir alias
  # its source dir and cause the patched CMakeLists.txt to overwrite the source
  # include/metis.h on every reconfigure (forcing a full METIS rebuild).
  add_subdirectory("${GKLIB_DIR}" "${PROJECT_BINARY_DIR}/GKlib")
  add_subdirectory("${METIS_DIR}" "${PROJECT_BINARY_DIR}/METIS")

  set(TPL_METIS_LIBRARIES "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/libmetis.a;${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/libGKlib.a" CACHE STRING "METIS static libraries" FORCE)
  set(TPL_METIS_INCLUDE_DIRS "${PROJECT_BINARY_DIR}/METIS/include" CACHE PATH "METIS include directory" FORCE)

  message(STATUS "Applying patches to SUPERLU...")
  file(GLOB_RECURSE SUPERLU_PATCH_FILES RELATIVE "${SUPERLU_PATCHES_DIR}" "${SUPERLU_PATCHES_DIR}/*")

  foreach(file IN LISTS SUPERLU_PATCH_FILES)
    set(src "${SUPERLU_PATCHES_DIR}/${file}")
    set(dst "${SUPERLU_DIR}/${file}")

    get_filename_component(dst_dir "${dst}" DIRECTORY)
    file(MAKE_DIRECTORY "${dst_dir}")

    message(STATUS "  Patching: ${file}")
    file(COPY "${src}" DESTINATION "${dst_dir}")
  endforeach()

  message(STATUS "SuperLU patched.")
endif()

set(enable_examples OFF FORCE)
set(enable_tests OFF FORCE)
set(enable_doc OFF FORCE)

# -----------------------------------------------------------------------
# BLAS / LAPACK provider resolution
#
# Single user-facing knob `MYSTRAN_BLAS_LAPACK` with three values:
# AUTO      (default) - try system BLAS+LAPACK; fall back to building
# Reference-LAPACK from the bundled submodule
# SYSTEM              - require system BLAS+LAPACK; FATAL_ERROR if missing
# EMBEDDED            - always build Reference-LAPACK + its reference BLAS
# from the Source/lapack submodule
#
# The mode is resolved here, BEFORE add_subdirectory() for SuperLU /
# SuperLU_MT, so we can forward the result to those subprojects via
# `TPL_BLAS_LIBRARIES` and `enable_internal_blaslib` and avoid duplicate
# symbols when statically linking on Windows (no symbol interposition).
#
# This replaces the older `MYSTRAN_BLAS` knob (which only governed BLAS,
# while LAPACK was always embedded under Source/Modules/LAPACK/*.f). The
# embedded LAPACK tree has been deleted; MYSTRAN now always links a real
# BLAS+LAPACK provider, and the only remaining MYSTRAN-specific routines
# (the handful of mathematically modified ones) live in
# Source/Modules/MYSTRAN_LAPACK_EXT/.
# -----------------------------------------------------------------------

# Backwards-compat: legacy `-Denable_internal_blaslib=...` and
# `-DMYSTRAN_BLAS=...` map to the new variable on first configure.
if(DEFINED enable_internal_blaslib AND NOT _MYSTRAN_BLAS_LAPACK_LEGACY_MAPPED)
  if(enable_internal_blaslib)
    set(_LEGACY_MODE "EMBEDDED")
  else()
    set(_LEGACY_MODE "SYSTEM")
  endif()

  message(DEPRECATION
    "'enable_internal_blaslib' is deprecated; "
    "use -DMYSTRAN_BLAS_LAPACK=${_LEGACY_MODE} instead.")
  set(MYSTRAN_BLAS_LAPACK "${_LEGACY_MODE}"
    CACHE STRING "BLAS+LAPACK provider: AUTO, SYSTEM, or EMBEDDED" FORCE)
  set(_MYSTRAN_BLAS_LAPACK_LEGACY_MAPPED TRUE CACHE INTERNAL "")
endif()

if(DEFINED MYSTRAN_BLAS AND NOT _MYSTRAN_BLAS_LAPACK_LEGACY_MAPPED)
  message(DEPRECATION
    "'MYSTRAN_BLAS' is deprecated; "
    "use -DMYSTRAN_BLAS_LAPACK=${MYSTRAN_BLAS} instead.")
  set(MYSTRAN_BLAS_LAPACK "${MYSTRAN_BLAS}"
    CACHE STRING "BLAS+LAPACK provider: AUTO, SYSTEM, or EMBEDDED" FORCE)
  set(_MYSTRAN_BLAS_LAPACK_LEGACY_MAPPED TRUE CACHE INTERNAL "")
endif()

set(MYSTRAN_BLAS_LAPACK "AUTO"
  CACHE STRING "BLAS+LAPACK provider: AUTO, SYSTEM, or EMBEDDED")
set_property(CACHE MYSTRAN_BLAS_LAPACK PROPERTY STRINGS AUTO SYSTEM EMBEDDED)

if(NOT MYSTRAN_BLAS_LAPACK MATCHES "^(AUTO|SYSTEM|EMBEDDED)$")
  message(FATAL_ERROR
    "MYSTRAN_BLAS_LAPACK must be AUTO, SYSTEM, or EMBEDDED "
    "(got '${MYSTRAN_BLAS_LAPACK}').")
endif()

set(_MYSTRAN_HAVE_SYSTEM_BLAS_LAPACK FALSE)

if(NOT MYSTRAN_BLAS_LAPACK STREQUAL "EMBEDDED")
  if(TPL_BLAS_LIBRARIES)
    set(BLAS_FOUND TRUE)
    set(BLAS_LIBRARIES "${TPL_BLAS_LIBRARIES}")
  else()
    find_package(BLAS QUIET)
  endif()

  if(BLAS_FOUND)
    find_package(LAPACK QUIET)
  endif()

  if(BLAS_FOUND AND LAPACK_FOUND)
    set(_MYSTRAN_HAVE_SYSTEM_BLAS_LAPACK TRUE)
  endif()
endif()

if(MYSTRAN_BLAS_LAPACK STREQUAL "SYSTEM" AND NOT _MYSTRAN_HAVE_SYSTEM_BLAS_LAPACK)
  message(FATAL_ERROR
    "MYSTRAN_BLAS_LAPACK=SYSTEM was requested but a system BLAS+LAPACK "
    "could not be located.\n"
    "Install OpenBLAS (Debian/Ubuntu: 'libopenblas-dev'; "
    "MSYS2/MinGW64: 'mingw-w64-x86_64-openblas') or set TPL_BLAS_LIBRARIES "
    "manually. Use -DMYSTRAN_BLAS_LAPACK=EMBEDDED to build the bundled "
    "Reference-LAPACK submodule instead.")
endif()

if(_MYSTRAN_HAVE_SYSTEM_BLAS_LAPACK)
  set(_MYSTRAN_BLAS_LAPACK_MODE "SYSTEM")
else()
  set(_MYSTRAN_BLAS_LAPACK_MODE "EMBEDDED")
endif()

# In EMBEDDED mode, build Reference-LAPACK (BLAS + LAPACK) from the
# submodule under Source/lapack. Configure it as a quiet, minimal,
# static-only sub-build before any add_subdirectory() that consumes it.
if(_MYSTRAN_BLAS_LAPACK_MODE STREQUAL "EMBEDDED")
  if(NOT EXISTS "${MYSTRAN_LAPACK_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
      "The Reference-LAPACK submodule was not downloaded. "
      "Run: git submodule update --init --recursive")
  endif()

  set(LAPACK_PATCHES_DIR "${SUBMODULES_DIR}/lapack_patches")
  message(STATUS "Applying patches to Reference-LAPACK...")
  file(GLOB_RECURSE LAPACK_PATCH_FILES RELATIVE "${LAPACK_PATCHES_DIR}" "${LAPACK_PATCHES_DIR}/*")

  foreach(file IN LISTS LAPACK_PATCH_FILES)
    set(src "${LAPACK_PATCHES_DIR}/${file}")
    set(dst "${MYSTRAN_LAPACK_DIR}/${file}")

    get_filename_component(dst_dir "${dst}" DIRECTORY)
    file(MAKE_DIRECTORY "${dst_dir}")

    message(STATUS "  Patching: ${file}")
    file(COPY "${src}" DESTINATION "${dst_dir}")
  endforeach()

  message(STATUS "Reference-LAPACK patched.")

  # Only build what MYSTRAN actually uses: double-precision real BLAS
  # and LAPACK. Everything else off.
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(BUILD_DEPRECATED OFF CACHE BOOL "" FORCE)
  set(BUILD_SINGLE OFF CACHE BOOL "" FORCE)
  set(BUILD_DOUBLE ON CACHE BOOL "" FORCE)
  set(BUILD_COMPLEX OFF CACHE BOOL "" FORCE)
  set(BUILD_COMPLEX16 OFF CACHE BOOL "" FORCE)
  set(LAPACKE OFF CACHE BOOL "" FORCE)
  set(CBLAS OFF CACHE BOOL "" FORCE)
  set(USE_OPTIMIZED_BLAS OFF CACHE BOOL "" FORCE)
  set(USE_OPTIMIZED_LAPACK OFF CACHE BOOL "" FORCE)

  # Do not let prior provider discovery bypass the embedded build on reconfigure.
  set(BLAS_LIBRARIES "")
  set(LAPACK_LIBRARIES "")
  set(BLAS_FOUND FALSE)
  set(LATESTLAPACK_FOUND FALSE)
  add_subdirectory("${MYSTRAN_LAPACK_DIR}" "${PROJECT_BINARY_DIR}/lapack" EXCLUDE_FROM_ALL)

  set(BLAS_LIBRARIES "blas")
  set(LAPACK_LIBRARIES "lapack")
  set(BLAS_LIBRARIES "blas" CACHE STRING "" FORCE)
  set(LAPACK_LIBRARIES "lapack" CACHE STRING "" FORCE)
endif()

# Forward the resolution to SuperLU / SuperLU_MT subprojects. They read
# `enable_internal_blaslib`, `TPL_ENABLE_INTERNAL_BLASLIB`, and
# `TPL_BLAS_LIBRARIES` to decide whether to build their bundled CBLAS.
if(_MYSTRAN_BLAS_LAPACK_MODE STREQUAL "SYSTEM")
  set(TPL_BLAS_LIBRARIES "${BLAS_LIBRARIES}"
    CACHE FILEPATH "Set from MYSTRAN BLAS resolution." FORCE)
  set(enable_internal_blaslib OFF
    CACHE BOOL "Use system BLAS for SuperLU" FORCE)
  set(TPL_ENABLE_INTERNAL_BLASLIB OFF
    CACHE BOOL "Use system BLAS for SuperLU" FORCE)
  message(STATUS
    "MYSTRAN BLAS+LAPACK mode: SYSTEM (BLAS=${BLAS_LIBRARIES}, LAPACK=${LAPACK_LIBRARIES})")
else()
  # Embedded mode: SuperLU links the Reference-LAPACK `blas` target
  # (Fortran BLAS only). SuperLU's bundled CBLAS shim is no longer
  # needed because Reference-LAPACK's blas archive provides every BLAS
  # symbol SuperLU references.
  set(enable_internal_blaslib OFF
    CACHE BOOL "Use embedded Reference-LAPACK BLAS for SuperLU" FORCE)
  set(TPL_ENABLE_INTERNAL_BLASLIB OFF
    CACHE BOOL "Use embedded Reference-LAPACK BLAS for SuperLU" FORCE)
  set(TPL_BLAS_LIBRARIES "blas"
    CACHE STRING "Set from MYSTRAN BLAS resolution." FORCE)
  message(STATUS
    "MYSTRAN BLAS+LAPACK mode: EMBEDDED (Reference-LAPACK submodule)")
endif()

if(USE_SUPERLU_MT)
  message(STATUS "Will link against faster SuperLU_MT.")
  set(PLAT "_OPENMP" CACHE STRING "threading flavor _PTHREAD/_OPENMP" FORCE)
  message(STATUS "Applying patches to SuperLU_MT...")
  file(GLOB_RECURSE PATCH_FILES RELATIVE "${SUPERLU_MT_PATCHES_DIR}" "${SUPERLU_MT_PATCHES_DIR}/*")

  foreach(file IN LISTS PATCH_FILES)
    set(src "${SUPERLU_MT_PATCHES_DIR}/${file}")
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

  if(TPL_ENABLE_METISLIB)
    if(PLAT STREQUAL "_PTHREAD")
      add_dependencies(superlu_mt_PTHREAD metis GKlib)
    elseif(PLAT STREQUAL "_OPENMP")
      add_dependencies(superlu_mt_OPENMP metis GKlib)
    endif()
  endif()

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

  if(TPL_ENABLE_METISLIB)
    add_dependencies(superlu metis GKlib)
  endif()

  set(SLU_DRIVER "${SUPERLU_DIR}/FORTRAN/c_fortran_dgssv.c")
endif()

# set some extra vars for MSYS builds to make the binary portable
if(WIN32)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
endif()

# collect all fortran source files
file(GLOB_RECURSE ALL_FORTRAN_FILES
  "${CMAKE_SOURCE_DIR}/*.f"
  "${CMAKE_SOURCE_DIR}/*.F"
  "${CMAKE_SOURCE_DIR}/*.f03"
  "${CMAKE_SOURCE_DIR}/*.F03"
)

# prepare the main executable, linked against the specifics and the m
# it appears utils used to be a module, but that is no longer the case?
# file(GLOB UTIL_FILES "${CMAKE_SOURCE_DIR}/UTIL/*.f*")
# file(GLOB MAIN_FILES "${CMAKE_SOURCE_DIR}/MAIN/*.[fF]*")
add_executable(
  mystran
  ${ALL_FORTRAN_FILES}
  ${SLU_DRIVER}
)

# determine which superlu variant to link against
if(USE_SUPERLU_MT)
  if(PLAT STREQUAL "_PTHREAD")
    message(WARNING "We recommend using OpenMP, not pthread!")
    target_link_libraries(mystran superlu_mt_PTHREAD)
  elseif(PLAT STREQUAL "_OPENMP")
    target_link_libraries(mystran superlu_mt_OPENMP)
  endif()

  target_compile_definitions(mystran PRIVATE USE_SUPERLU_MT)

  if(TPL_ENABLE_METISLIB)
    target_compile_definitions(mystran PRIVATE HAVE_METIS)
  endif()
else()
  target_link_libraries(mystran superlu)
endif()

# Link LAPACK first, then BLAS: LAPACK depends on BLAS so the linker
# must see BLAS symbols last when resolving LAPACK references.
target_link_libraries(mystran ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

set_target_properties(
  mystran PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

if(WIN32)
  target_sources(mystran PRIVATE mystran.rc)
endif()

# -----------------------------------------------------------------------
# Static library tracking
#
# For every library that is statically linked into mystran we:
# 1. Record a compile-time preprocessor symbol (_STATIC_<NAME>)
# 2. Record a human-readable filename for the status message
#
# After populating the lists we apply the definitions, print a summary,
# and warn for any symbol not referenced in STARTUP_SPECIALS.f03.
# -----------------------------------------------------------------------
set(_MYSTRAN_STATIC_DEFS "")
set(_MYSTRAN_STATIC_NAMES "")

if(USE_SUPERLU_MT)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_SUPERLU_MT)

  if(WIN32)
    list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_WINPTHREAD)
    list(APPEND _MYSTRAN_STATIC_NAMES libwinpthread.a)
  else()
    list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_PTHREAD)

    # no name because it's part of glibc
  endif()

  if(PLAT STREQUAL "_OPENMP")
    list(APPEND _MYSTRAN_STATIC_NAMES libsuperlu_mt_OPENMP.a)

    if(WIN32)
      # On Windows the gomp runtime is force-linked as a static archive
      list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_GOMP)
      list(APPEND _MYSTRAN_STATIC_NAMES libgomp.a)
    endif()
  elseif(PLAT STREQUAL "_PTHREAD")
    list(APPEND _MYSTRAN_STATIC_NAMES libsuperlu_mt_PTHREAD.a)
  endif()
else()
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_SUPERLU)
  list(APPEND _MYSTRAN_STATIC_NAMES libsuperlu.a)
endif()

if(TPL_ENABLE_METISLIB)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_METIS)
  list(APPEND _MYSTRAN_STATIC_NAMES libmetis.a)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_GKLIB)
  list(APPEND _MYSTRAN_STATIC_NAMES libGKlib.a)

  if(NOT WIN32)
    list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_REGEX)

    # only on linux because on windows we build embedded gkegex.
    # no name because it's part of glibc.
  endif()
endif()

# BLAS/LAPACK provider tracking.
#
# EMBEDDED: the Reference-LAPACK submodule provides static libblas.a
# and liblapack.a archives, and its own LICENSE is bundled.
# SYSTEM:   we statically link a system BLAS+LAPACK. We try to detect
# the vendor (OpenBLAS, MKL, BLIS, ATLAS, Apple Accelerate,
# Arm Performance Libraries) from BLAS_LIBRARIES/LAPACK_LIBRARIES
# so the corresponding license text from dependencies_licenses/ can
# be embedded in the binary. If detection fails OR the license
# file is missing, we WARN the user that the resulting binary
# may not be legally redistributable.
set(_SYS_BLAS_TAG "")
set(_SYS_BLAS_LIBNAME "")
set(_SYS_BLAS_DESC "")
set(_SYS_BLAS_LICENSE_PATH "")

if(_MYSTRAN_BLAS_LAPACK_MODE STREQUAL "EMBEDDED")
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_BLAS)
  list(APPEND _MYSTRAN_STATIC_NAMES libblas.a)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_LAPACK)
  list(APPEND _MYSTRAN_STATIC_NAMES liblapack.a)
else()
  # Use the concatenated link-line strings as the signature to match.
  string(TOLOWER "${BLAS_LIBRARIES};${LAPACK_LIBRARIES}" _SYS_BLAS_SIG)

  if(_SYS_BLAS_SIG MATCHES "openblas")
    set(_SYS_BLAS_TAG "OPENBLAS")
    set(_SYS_BLAS_LIBNAME "OpenBLAS")
    set(_SYS_BLAS_DESC "system BLAS+LAPACK")
    set(_SYS_BLAS_LICENSE_PATH "dependencies_licenses/openblas.txt")
  elseif(_SYS_BLAS_SIG MATCHES "mkl")
    set(_SYS_BLAS_TAG "MKL")
    set(_SYS_BLAS_LIBNAME "Intel MKL")
    set(_SYS_BLAS_DESC "system BLAS+LAPACK")
    set(_SYS_BLAS_LICENSE_PATH "dependencies_licenses/mkl.txt")
  elseif(_SYS_BLAS_SIG MATCHES "blis")
    set(_SYS_BLAS_TAG "BLIS")
    set(_SYS_BLAS_LIBNAME "BLIS")
    set(_SYS_BLAS_DESC "system BLAS")
    set(_SYS_BLAS_LICENSE_PATH "dependencies_licenses/blis.txt")
  elseif(_SYS_BLAS_SIG MATCHES "atlas")
    set(_SYS_BLAS_TAG "ATLAS")
    set(_SYS_BLAS_LIBNAME "ATLAS")
    set(_SYS_BLAS_DESC "system BLAS+LAPACK")
    set(_SYS_BLAS_LICENSE_PATH "dependencies_licenses/atlas.txt")
  elseif(_SYS_BLAS_SIG MATCHES "blas")
    # Generic Reference-BLAS / netlib install on the system.
    set(_SYS_BLAS_TAG "REFBLAS")
    set(_SYS_BLAS_LIBNAME "Netlib LAPACK")
    set(_SYS_BLAS_DESC "system BLAS+LAPACK")
    set(_SYS_BLAS_LICENSE_PATH "dependencies_licenses/netlib-blas-lapack.txt")
  endif()

  if(_SYS_BLAS_TAG STREQUAL "")
    message(WARNING
      "Could not identify the system BLAS/LAPACK vendor from "
      "BLAS_LIBRARIES='${BLAS_LIBRARIES}', LAPACK_LIBRARIES='${LAPACK_LIBRARIES}'.\n"
      "  No license text will be embedded for the BLAS/LAPACK provider, "
      "which may make the resulting binary legally undistributable.\n"
      "  Add a matching entry to the system-BLAS detection block in "
      "CMakeLists.txt and drop the license file under dependencies_licenses/.")
  else()
    list(APPEND _MYSTRAN_STATIC_DEFS "_STATIC_${_SYS_BLAS_TAG}")
    list(APPEND _MYSTRAN_STATIC_NAMES "${_SYS_BLAS_LIBNAME}")

    # Provider-specific runtime hooks (e.g. thread-count limiter in
    # Source/MAIN/SET_BLAS_THREADS.F03). Defined here so it's only
    # active when the corresponding BLAS is actually linked.
    if(_SYS_BLAS_TAG STREQUAL "OPENBLAS")
      target_compile_definitions(mystran PRIVATE MYSTRAN_USE_OPENBLAS)
    endif()

    if(NOT EXISTS "${PROJECT_SOURCE_DIR}/${_SYS_BLAS_LICENSE_PATH}")
      message(WARNING
        "System BLAS/LAPACK detected as ${_SYS_BLAS_LIBNAME}, but the "
        "expected license file is missing:\n"
        "  ${PROJECT_SOURCE_DIR}/${_SYS_BLAS_LICENSE_PATH}\n"
        "  The binary will be built without an embedded license for "
        "${_SYS_BLAS_LIBNAME}, which may make it legally undistributable.")
    endif()
  endif()
endif()

# On Windows the -static linker flag pulls in Fortran and C runtimes too
if(WIN32)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_GFORTRAN)
  list(APPEND _MYSTRAN_STATIC_NAMES libgfortran.a)
  list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_GCC)
  list(APPEND _MYSTRAN_STATIC_NAMES libgcc.a)
endif()

# Apply preprocessor definitions to the mystran compilation
foreach(_def IN LISTS _MYSTRAN_STATIC_DEFS)
  target_compile_definitions(mystran PRIVATE ${_def})
endforeach()

# Print a human-readable summary at configure time
list(JOIN _MYSTRAN_STATIC_NAMES ", " _STATIC_NAMES_STR)
message(STATUS "Static libraries linked in: ${_STATIC_NAMES_STR}")

# ARPACK is always embedded under Source/Modules/ARPACK.
list(APPEND _MYSTRAN_STATIC_DEFS _STATIC_ARPACK)

# -----------------------------------------------------------------------
# License file associations
#
# Map each _STATIC_<NAME> suffix to a license file path (relative to
# PROJECT_SOURCE_DIR).  CMake will:
# - define  _STATIC_<NAME>_LICENSE  for the Fortran preprocessor
# - generate a SUBROUTINE PRINT_LICENSE_<NAME>(IUNIT) that writes the
# license text line-by-line to whatever unit number is passed in
#
# Add new entries here whenever a new statically-linked library is added.
# -----------------------------------------------------------------------
set(_MYSTRAN_LICENSE_MAP
  "SUPERLU:submodules/superlu/License.txt:SuperLU:sparse solver"
  "SUPERLU_MT:submodules/superlu_mt/License.txt:SuperLU_MT:multi-threaded sparse solver"
  "METIS:submodules/metis/LICENSE:METIS:graph partitioner for SuperLU"
  "GKLIB:submodules/GKlib/LICENSE.txt:GKlib:utility library for METIS"
  "ARPACK:dependencies_licenses/arpack.txt:ARPACK:Lanczos eigensolver"
  "GOMP:dependencies_licenses/gpl-3.0-with-gcc-exception.txt:libgomp:GNU OpenMP runtime"
  "GFORTRAN:dependencies_licenses/gpl-3.0-with-gcc-exception.txt:libgfortran:GNU Fortran runtime"
  "GCC:dependencies_licenses/gpl-3.0-with-gcc-exception.txt:libgcc:GNU C runtime"
  "WINPTHREAD:dependencies_licenses/libwinpthread.txt:libwinpthread:Windows POSIX threads, used by SuperLU_MT"
  "PTHREAD:dependencies_licenses/lgpl-2.1.txt:libpthread:POSIX threads, used by SuperLU_MT"
  "REGEX:dependencies_licenses/lgpl-2.1.txt:libregex:POSIX regexes, used by GKlib"
)

# BLAS/LAPACK license entry depends on which provider was selected.
# EMBEDDED -> Reference-LAPACK submodule ships its own LICENSE.
# SYSTEM   -> use the detected vendor's license file, if known.
if(_MYSTRAN_BLAS_LAPACK_MODE STREQUAL "EMBEDDED")
  list(APPEND _MYSTRAN_LICENSE_MAP
    "BLAS:submodules/lapack/LICENSE:BLAS:Reference-LAPACK submodule (embedded BLAS)"
    "LAPACK:submodules/lapack/LICENSE:LAPACK:Reference-LAPACK submodule (embedded LAPACK)"
  )
elseif(NOT _SYS_BLAS_TAG STREQUAL "")
  list(APPEND _MYSTRAN_LICENSE_MAP
    "${_SYS_BLAS_TAG}:${_SYS_BLAS_LICENSE_PATH}:${_SYS_BLAS_LIBNAME}:${_SYS_BLAS_DESC}"
  )
endif()

# Generate STATIC_LICENSES_GENERATED.F03 in the build directory
set(_GEN_LIC_FILE "${CMAKE_SOURCE_DIR}/INCLUDE/STATIC_LICENSES_GENERATED.F03.inc")
set(_GEN_LIC_CONTENT "! AUTO-GENERATED by CMakeLists.txt -- DO NOT EDIT\n")
set(_MYSTRAN_LICENSED_NAMES "")

# ---- MYSTRAN's own license (always emitted, no #ifdef guard) ----
set(_MYSTRAN_OWN_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE.txt")

if(EXISTS "${_MYSTRAN_OWN_LICENSE}")
  file(STRINGS "${_MYSTRAN_OWN_LICENSE}" _mlic_lines NEWLINE_CONSUME)
  string(REPLACE "\n" ";" _mlic_lines "${_mlic_lines}")
  string(APPEND _GEN_LIC_CONTENT
    "\nSUBROUTINE PRINT_LICENSE_MYSTRAN(IUNIT)\n"
    "  IMPLICIT NONE\n"
    "  INTEGER, INTENT(IN) :: IUNIT\n"
  )

  foreach(_line IN LISTS _mlic_lines)
    string(REPLACE "'" "''" _line_esc "${_line}")
    string(REPLACE "\\" "\\\\" _line_esc "${_line_esc}")
    string(APPEND _GEN_LIC_CONTENT
      "  WRITE(IUNIT,'(A)') '${_line_esc}'\n"
    )
  endforeach()

  string(APPEND _GEN_LIC_CONTENT
    "END SUBROUTINE PRINT_LICENSE_MYSTRAN\n"
  )
else()
  message(WARNING "MYSTRAN license file not found: ${_MYSTRAN_OWN_LICENSE}")
endif()

foreach(_entry IN LISTS _MYSTRAN_LICENSE_MAP)
  string(REPLACE ":" ";" _pair "${_entry}")
  list(GET _pair 0 _lname)
  list(GET _pair 1 _lpath)
  list(LENGTH _pair _pair_len)

  if(_pair_len GREATER 2)
    list(GET _pair 2 _llibname)
  else()
    set(_llibname "${_lname}")
  endif()

  if(_pair_len GREATER 3)
    list(GET _pair 3 _ldesc)
  else()
    set(_ldesc "")
  endif()

  set(_MYSTRAN_LIBNAME_${_lname} "${_llibname}")
  set(_MYSTRAN_DESC_${_lname} "${_ldesc}")
  set(_lfile "${PROJECT_SOURCE_DIR}/${_lpath}")

  if(EXISTS "${_lfile}")
    # Read the license file line by line
    file(STRINGS "${_lfile}" _lic_lines NEWLINE_CONSUME)

    # file(STRINGS) may return the whole file as one item when using
    # NEWLINE_CONSUME; split on actual newlines ourselves
    string(REPLACE "\n" ";" _lic_lines "${_lic_lines}")

    string(APPEND _GEN_LIC_CONTENT
      "\nSUBROUTINE PRINT_LICENSE_${_lname}(IUNIT)\n"
      "  IMPLICIT NONE\n"
      "  INTEGER, INTENT(IN) :: IUNIT\n"
    )

    foreach(_line IN LISTS _lic_lines)
      # Escape Fortran string delimiters: ' -> ''
      string(REPLACE "'" "''" _line_esc "${_line}")

      # Escape backslashes that CMake would otherwise swallow
      string(REPLACE "\\" "\\\\" _line_esc "${_line_esc}")
      string(APPEND _GEN_LIC_CONTENT
        "  WRITE(IUNIT,'(A)') '${_line_esc}'\n"
      )
    endforeach()

    string(APPEND _GEN_LIC_CONTENT
      "END SUBROUTINE PRINT_LICENSE_${_lname}\n"
    )

    list(APPEND _MYSTRAN_LICENSED_NAMES "${_lname}")
  else()
    message(WARNING
      "License file for _STATIC_${_lname} not found: ${_lfile}\n"
      "  '_STATIC_${_lname}_LICENSE' will NOT be defined."
    )
  endif()
endforeach()

# Append PRINT_ALL_LICENSES: MYSTRAN first, then delegates to each PRINT_LICENSE_* under its #ifdef guard
string(APPEND _GEN_LIC_CONTENT
  "\nSUBROUTINE PRINT_ALL_LICENSES(IUNIT)\n"
  "  IMPLICIT NONE\n"
  "  INTEGER, INTENT(IN) :: IUNIT\n"

  # "  WRITE(IUNIT,'(A)') '========================================================'\n"
  # "  WRITE(IUNIT,'(A)') ' MYSTRAN '\n"
  "  WRITE(IUNIT,'(A)') '========================================================'\n"
  "  CALL PRINT_LICENSE_MYSTRAN(IUNIT)\n"
  "  WRITE(IUNIT,'(A)') ' '\n"
)

foreach(_lname IN LISTS _MYSTRAN_LICENSED_NAMES)
  set(_llibname "${_MYSTRAN_LIBNAME_${_lname}}")
  set(_ldesc "${_MYSTRAN_DESC_${_lname}}")
  string(APPEND _GEN_LIC_CONTENT
    "#ifdef _STATIC_${_lname}\n"
    "  WRITE(IUNIT,'(A)') '========================================================'\n"
    "  WRITE(IUNIT,'(A)') ' ${_llibname} -- ${_ldesc}'\n"
    "  WRITE(IUNIT,'(A)') '========================================================'\n"
    "  CALL PRINT_LICENSE_${_lname}(IUNIT)\n"
    "  WRITE(IUNIT,'(A)') ' '\n"
    "#endif\n"
  )
endforeach()

string(APPEND _GEN_LIC_CONTENT
  "END SUBROUTINE PRINT_ALL_LICENSES\n"
)

file(WRITE "${_GEN_LIC_FILE}" "${_GEN_LIC_CONTENT}")
message(STATUS "Generated license subroutines: ${_GEN_LIC_FILE}")

# Warn if any statically-linked library lacks a license mapping.
# Runtime and utility libs with no separate OSS license file are exempt.
set(_LICENSE_EXEMPT_LIBS)

foreach(_def IN LISTS _MYSTRAN_STATIC_DEFS)
  string(REPLACE "_STATIC_" "" _lname "${_def}")

  if(NOT _lname IN_LIST _MYSTRAN_LICENSED_NAMES)
    if(NOT _lname IN_LIST _LICENSE_EXEMPT_LIBS)
      message(WARNING
        "No license file is mapped for '${_def}'. "
        "Add an entry to _MYSTRAN_LICENSE_MAP in CMakeLists.txt."
      )
    endif()
  endif()
endforeach()

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_Fortran_FLAGS "-Wall -Wextra \
    -Wno-unused-variable -Wno-unused-label -Wno-unused-parameter -Wno-tabs \
    -Wno-compare-reals -Wno-character-truncation -Wno-unused-dummy-argument \
    -Wmaybe-uninitialized -Wrealloc-lhs -fbounds-check -fcheck=all \
    -ffree-line-length-none -fno-diagnostics-show-caret")
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

# -----------------------------------------------------------------------
# Build metadata generated file
#
# Generates BUILD_INFO_GENERATED.F03.inc with:
# PRINT_BUILD_CONSTANTS(IUNIT) - timestamp, compiler, flags, solver info
# PRINT_STATIC_LIB_LIST(IUNIT) - names of all statically linked libraries
# -----------------------------------------------------------------------

# Build timestamp captured at configure time (local clock)
string(TIMESTAMP _BUILD_TIMESTAMP "%Y-%m-%d %H:%M:%S UTC" UTC)

set(_BUILD_TYPE_DISPLAY "${CMAKE_BUILD_TYPE}")

if(NOT _BUILD_TYPE_DISPLAY)
  set(_BUILD_TYPE_DISPLAY "(not specified)")
endif()

# Effective Fortran flags = base flags + build-type-specific flags
string(TOUPPER "${CMAKE_BUILD_TYPE}" _BT_UPPER)
set(_FC_FLAGS "${CMAKE_Fortran_FLAGS}")

if(_BT_UPPER AND DEFINED CMAKE_Fortran_FLAGS_${_BT_UPPER})
  string(APPEND _FC_FLAGS " ${CMAKE_Fortran_FLAGS_${_BT_UPPER}}")
endif()

string(JOIN " " _MYSTRAN_CPU_FLAGS_DISPLAY ${_MYSTRAN_CPU_FLAGS})
string(APPEND _FC_FLAGS " ${_MYSTRAN_CPU_FLAGS_DISPLAY}")

string(REPLACE "\n" " " _FC_FLAGS "${_FC_FLAGS}")
string(REGEX REPLACE " +" " " _FC_FLAGS "${_FC_FLAGS}")
string(STRIP "${_FC_FLAGS}" _FC_FLAGS)
string(REPLACE "'" "''" _FC_FLAGS "${_FC_FLAGS}")

# Effective C flags
set(_CC_FLAGS "${CMAKE_C_FLAGS}")

if(_BT_UPPER AND DEFINED CMAKE_C_FLAGS_${_BT_UPPER})
  string(APPEND _CC_FLAGS " ${CMAKE_C_FLAGS_${_BT_UPPER}}")
endif()

string(APPEND _CC_FLAGS " ${_MYSTRAN_CPU_FLAGS_DISPLAY}")

string(REPLACE "\n" " " _CC_FLAGS "${_CC_FLAGS}")
string(REGEX REPLACE " +" " " _CC_FLAGS "${_CC_FLAGS}")
string(STRIP "${_CC_FLAGS}" _CC_FLAGS)
string(REPLACE "'" "''" _CC_FLAGS "${_CC_FLAGS}")

# Compiler identity strings
set(_FC_DISPLAY "${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
set(_CC_DISPLAY "${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
string(REPLACE "'" "''" _FC_PATH "${CMAKE_Fortran_COMPILER}")
string(REPLACE "'" "''" _CC_PATH "${CMAKE_C_COMPILER}")

# Sparse solver variant
if(USE_SUPERLU_MT)
  if(PLAT STREQUAL "_OPENMP")
    set(_SLU_VARIANT "SuperLU_MT (OpenMP)")
  else()
    set(_SLU_VARIANT "SuperLU_MT (pthreads)")
  endif()
else()
  set(_SLU_VARIANT "SuperLU (single-threaded)")
endif()

# BLAS / LAPACK source
if(_MYSTRAN_BLAS_LAPACK_MODE STREQUAL "SYSTEM")
  set(_BLAS_INFO "System: ${BLAS_LIBRARIES}")
  set(_LAPACK_INFO "System: ${LAPACK_LIBRARIES}")
else()
  set(_BLAS_INFO "Embedded Reference-LAPACK submodule (libblas.a)")
  set(_LAPACK_INFO "Embedded Reference-LAPACK submodule (liblapack.a)")
endif()

string(REPLACE "'" "''" _BLAS_INFO "${_BLAS_INFO}")
string(REPLACE "'" "''" _LAPACK_INFO "${_LAPACK_INFO}")

# METIS
if(TPL_ENABLE_METISLIB)
  set(_METIS_INFO "Enabled")
else()
  set(_METIS_INFO "Disabled")
endif()

set(_GEN_BLD_FILE "${CMAKE_SOURCE_DIR}/INCLUDE/BUILD_INFO_GENERATED.F03.inc")
set(_GEN_BLD_CONTENT "! AUTO-GENERATED by CMakeLists.txt -- DO NOT EDIT\n")

# PRINT_BUILD_CONSTANTS subroutine
string(APPEND _GEN_BLD_CONTENT
  "\nSUBROUTINE PRINT_BUILD_CONSTANTS(IUNIT)\n"
  "  IMPLICIT NONE\n"
  "  INTEGER, INTENT(IN) :: IUNIT\n"
  "  WRITE(IUNIT,'(A)') ' '\n"
  "  WRITE(IUNIT,'(A)') ' Build information:'\n"
  "  WRITE(IUNIT,'(A)') ' ------------------'\n"
  "  WRITE(IUNIT,'(A)') '  Build date (UTC)    : ${_BUILD_TIMESTAMP}'\n"
  "  WRITE(IUNIT,'(A)') '  CMake build type    : ${_BUILD_TYPE_DISPLAY}'\n"
  "  WRITE(IUNIT,'(A)') ' '\n"
  "  WRITE(IUNIT,'(A)') '  Fortran compiler    : ${_FC_DISPLAY}'\n"
  "  WRITE(IUNIT,'(A)') '  Fortran path        : ${_FC_PATH}'\n"
  "  WRITE(IUNIT,'(A)') '  Fortran flags       : ${_FC_FLAGS}'\n"
  "  WRITE(IUNIT,'(A)') ' '\n"
  "  WRITE(IUNIT,'(A)') '  C compiler          : ${_CC_DISPLAY}'\n"
  "  WRITE(IUNIT,'(A)') '  C compiler path     : ${_CC_PATH}'\n"
  "  WRITE(IUNIT,'(A)') '  C flags             : ${_CC_FLAGS}'\n"
  "  WRITE(IUNIT,'(A)') ' '\n"
  "  WRITE(IUNIT,'(A)') '  Sparse solver       : ${_SLU_VARIANT}'\n"
  "  WRITE(IUNIT,'(A)') '  BLAS                : ${_BLAS_INFO}'\n"
  "  WRITE(IUNIT,'(A)') '  LAPACK              : ${_LAPACK_INFO}'\n"
  "  WRITE(IUNIT,'(A)') '  METIS partitioner   : ${_METIS_INFO}'\n"
  "END SUBROUTINE PRINT_BUILD_CONSTANTS\n"
)

# PRINT_STATIC_LIB_LIST subroutine
string(APPEND _GEN_BLD_CONTENT
  "\nSUBROUTINE PRINT_STATIC_LIB_LIST(IUNIT)\n"
  "  IMPLICIT NONE\n"
  "  INTEGER, INTENT(IN) :: IUNIT\n"
  "  WRITE(IUNIT,'(A)') ' Static libraries compiled into this MYSTRAN binary:'\n"
  "  WRITE(IUNIT,'(A)') ' ----------------------------------------------------'\n"
)

foreach(_def IN LISTS _MYSTRAN_STATIC_DEFS)
  string(REPLACE "_STATIC_" "" _suffix "${_def}")
  set(_libname "${_MYSTRAN_LIBNAME_${_suffix}}")
  set(_desc "${_MYSTRAN_DESC_${_suffix}}")

  if(NOT _libname)
    set(_libname "${_def}")
    set(_desc "(no description)")
  endif()

  string(REPLACE "'" "''" _libname_esc "${_libname}")
  string(REPLACE "'" "''" _desc_esc "${_desc}")

  string(APPEND _GEN_BLD_CONTENT
    "  WRITE(IUNIT,'(3X,A,T20,A)') '${_libname_esc}', '(${_desc_esc})'\n"
  )
endforeach()

string(APPEND _GEN_BLD_CONTENT
  "  WRITE(IUNIT,'(A)') ' '\n"
  "END SUBROUTINE PRINT_STATIC_LIB_LIST\n"
)

file(WRITE "${_GEN_BLD_FILE}" "${_GEN_BLD_CONTENT}")
message(STATUS "Generated build-info subroutines: ${_GEN_BLD_FILE}")

# issue a couple messages about compilation
include(ProcessorCount)
ProcessorCount(NCPU)
message(STATUS "You can now compile MYSTRAN with `cmake --build .`")

if(NOT NCPU EQUAL 0)
  message(STATUS "Compile faster by passing the -j${NCPU} flag.")
endif()
