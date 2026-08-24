# Standalone build for the multi-membrane test driver.
#
# Links libmem3dg's C++20 API directly (not the plain-C ABI), so it needs the
# same standard and the same geometry-central / Eigen headers the library was
# built against.  Point MEM3DG_BUILD at the directory used in BUILD.md step 1.
#
#   cmake -S tests/test-multi-membrane -B build-multi \
#     -DMEM3DG_ROOT=/path/to/Mem3DG -DMEM3DG_BUILD=/path/to/mem3dg-build
#   cmake --build build-multi -j8

cmake_minimum_required(VERSION 3.14)
project(multi_membrane_driver CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT DEFINED MEM3DG_ROOT)
  message(FATAL_ERROR "Set -DMEM3DG_ROOT=/abs/path/to/Mem3DG")
endif()
if(NOT DEFINED MEM3DG_BUILD)
  message(FATAL_ERROR "Set -DMEM3DG_BUILD=/abs/path/to/the/mem3dg/build/dir")
endif()

find_library(MEM3DG_LIBRARY NAMES mem3dg PATHS "${MEM3DG_BUILD}/lib" REQUIRED)

# geometry-central is VENDORED in the repo (libraries/), while Eigen is fetched
# by FetchContent into the libmem3dg build tree.  Two different places, hence
# two different variables.
set(_gc_src "${MEM3DG_ROOT}/libraries/geometry-central")
if(NOT EXISTS "${_gc_src}/include/geometrycentral")
  message(FATAL_ERROR "geometry-central headers not at ${_gc_src}")
endif()
set(_eigen_src "${MEM3DG_BUILD}/_deps/eigen-src")
if(NOT EXISTS "${_eigen_src}/Eigen")
  # FETCHCONTENT_BASE_DIR may point elsewhere; let the caller override
  message(WARNING "Eigen not at ${_eigen_src}; pass -DEIGEN_SRC=...")
endif()
if(DEFINED EIGEN_SRC)
  set(_eigen_src "${EIGEN_SRC}")
endif()

add_executable(multi_membrane_driver multi_membrane_driver.cc)

target_include_directories(multi_membrane_driver PRIVATE
    "${MEM3DG_ROOT}/include"
    "${MEM3DG_BUILD}/include"
    "${_gc_src}/include"
    "${_gc_src}/deps/nanort/include"
    "${_gc_src}/deps/nanoflann"
    "${_gc_src}/deps/happly"
    "${MEM3DG_ROOT}/libraries/pcg/include"
    "${MEM3DG_ROOT}/libraries/libigl/include"
    "${_eigen_src}"
)

# ---------------------------------------------------------------------------
# MEM3DG_WITH_NETCDF must MATCH the library, or the build is silently broken.
#
# That define adds members to System and Integrator (integrator.h:69 declares
# `MutableTrajFile mutableTrajFile' under it, system.h has three more sites).
# Compiling this target without it while libmem3dg has it gives the two
# translation units DIFFERENT struct layouts for the same class: field offsets
# shift, and code here reads the wrong bytes.  The failure is not a link error,
# it is corruption -- the observed symptom was ~MutableTrajFile() calling
# NcFile::close() on a pointer whose bit pattern was a double, crashing at exit
# AFTER the physics had printed plausible-looking numbers.
#
# So do not hardcode it: read it out of the library's own CMakeCache, which is
# the only authority on how libmem3dg was actually compiled.
# ---------------------------------------------------------------------------
set(_m3dg_cache "${MEM3DG_BUILD}/CMakeCache.txt")
if(EXISTS "${_m3dg_cache}")
  file(STRINGS "${_m3dg_cache}" _netcdf_line REGEX "^WITH_NETCDF:BOOL=")
  if(_netcdf_line MATCHES "=ON$")
    set(_m3dg_netcdf ON)
  else()
    set(_m3dg_netcdf OFF)
  endif()
  message(STATUS "libmem3dg was built with WITH_NETCDF=${_m3dg_netcdf}; matching it")
else()
  message(WARNING "no CMakeCache.txt in ${MEM3DG_BUILD}; assuming WITH_NETCDF=ON. "
                  "If that is wrong you will get struct-layout corruption, not a link error.")
  set(_m3dg_netcdf ON)
endif()

if(_m3dg_netcdf)
  target_compile_definitions(multi_membrane_driver PRIVATE MEM3DG_WITH_NETCDF)
  # the trajectory-file code is header-inlined here, so this target needs the
  # netcdf headers and library itself
  find_path(NETCDF_CXX4_INCLUDE_DIR ncFile.h
            HINTS /opt/homebrew/include /usr/local/include)
  find_library(NETCDF_CXX4_LIBRARY NAMES netcdf-cxx4
               HINTS /opt/homebrew/lib /usr/local/lib REQUIRED)
  target_include_directories(multi_membrane_driver PRIVATE "${NETCDF_CXX4_INCLUDE_DIR}")
endif()

# libmem3dg links geometry-central statically and does NOT re-export its
# symbols, while much of geometry-central is header-inlined into our own object
# file.  So the archive has to be linked here as well.
find_library(GC_LIBRARY NAMES geometry-central PATHS "${MEM3DG_BUILD}/lib" REQUIRED)

target_link_libraries(multi_membrane_driver PRIVATE
    "${MEM3DG_LIBRARY}" "${GC_LIBRARY}")

if(_m3dg_netcdf)
  target_link_libraries(multi_membrane_driver PRIVATE "${NETCDF_CXX4_LIBRARY}")
endif()

set_target_properties(multi_membrane_driver PROPERTIES
    BUILD_RPATH   "${MEM3DG_BUILD}/lib"
    INSTALL_RPATH "${MEM3DG_BUILD}/lib")
