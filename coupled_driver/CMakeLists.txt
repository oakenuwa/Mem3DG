# =============================================================================
# Standalone external driver coupling Cytosim and Mem3DG.
#
# Self-contained on purpose: it does not add anything to Cytosim's build, and
# it does not need Mem3DG's BUILD_COUPLED machinery.  It links
#   - Cytosim's static libraries, built separately
#   - libmem3dg.{so,dylib}, via its plain-C ABI only
#
# Because the only Mem3DG surface used is <mem3dg/capi/mem3dg_c.h>, this target
# compiles as C++11 -- the same standard Cytosim uses.  No C++20 anywhere.
#
# Configure with:
#   cmake -S coupled_driver -B build-driver \
#     -DCYTOSIM_DIR=/path/to/cytosim-lee-lab-2 \
#     -DCYTOSIM_BUILD=/path/to/cytosim/build   \
#     -DMEM3DG_ROOT=/path/to/Mem3DG            \
#     -DMEM3DG_LIB=/path/to/libmem3dg.dylib
# =============================================================================

cmake_minimum_required(VERSION 3.5)
project(coupled_driver CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- inputs ------------------------------------------------------------------
if(NOT DEFINED CYTOSIM_DIR)
  message(FATAL_ERROR "Set -DCYTOSIM_DIR=/abs/path/to/cytosim source root")
endif()
if(NOT DEFINED CYTOSIM_BUILD)
  message(FATAL_ERROR "Set -DCYTOSIM_BUILD=/abs/path/to/cytosim build dir")
endif()
if(NOT DEFINED MEM3DG_ROOT)
  message(FATAL_ERROR "Set -DMEM3DG_ROOT=/abs/path/to/Mem3DG")
endif()
if(NOT DEFINED MEM3DG_LIB)
  set(MEM3DG_LIB "${MEM3DG_ROOT}/build/lib/libmem3dg${CMAKE_SHARED_LIBRARY_SUFFIX}")
endif()
if(NOT EXISTS "${MEM3DG_LIB}")
  message(FATAL_ERROR "libmem3dg not found at ${MEM3DG_LIB}")
endif()

# Cytosim bakes the spatial dimension in at compile time; the driver's objects
# must agree with the libraries it links.  Mem3DG meshes are 3D.
set(DIMENSION 3 CACHE STRING "spatial dimension; must match the cytosim build")
add_compile_definitions(DIM=${DIMENSION})

# SpaceMem3DG is declared only under `DIM == 3 && HAS_MEM3DG' (see
# cytosim src/sim/spaces/space_mem3dg.h).  Cytosim's own build sets HAS_MEM3DG
# from -DMEM3DG_ROOT, but this driver locates the cytosim ARCHIVES rather than
# linking the CMake target, so that PUBLIC compile definition does not propagate
# and has to be repeated.  It is unconditional here: MEM3DG_ROOT and libmem3dg
# are both hard requirements above, so the condition always holds.  Without it
# coupled_driver.cc fails to compile on `SpaceMem3DG::setHandle'.
add_compile_definitions(HAS_MEM3DG=1)

# --- Cytosim static libraries ------------------------------------------------
# Built separately; we locate the archives rather than rebuilding them here.
find_library(CYTO_SIM_LIB   NAMES cytosim  PATHS "${CYTOSIM_BUILD}/tar" REQUIRED)
find_library(CYTO_MATH_LIB  NAMES cytomath PATHS "${CYTOSIM_BUILD}/tar" REQUIRED)
find_library(CYTO_BASE_LIB  NAMES cytobase PATHS "${CYTOSIM_BUILD}/tar" REQUIRED)

# Cytosim's headers are spread across these directories (see its src/sim/CMakeLists)
set(CYTO_INCLUDES math base sim sim/fibers sim/hands sim/singles
                  sim/couples sim/organizers sim/spaces)
list(TRANSFORM CYTO_INCLUDES PREPEND "${CYTOSIM_DIR}/src/")

# --- dependencies Cytosim's static libs leave unresolved ---------------------
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(Threads REQUIRED)
find_package(Curses REQUIRED)

# --- target ------------------------------------------------------------------
add_executable(coupled_driver coupled_driver.cc)

target_include_directories(coupled_driver PRIVATE
    "${CYTO_INCLUDES}"
    "${MEM3DG_ROOT}/include"      # for <mem3dg/capi/mem3dg_c.h>
    "${CURSES_INCLUDE_DIR}"
)

target_link_libraries(coupled_driver PRIVATE
    "${CYTO_SIM_LIB}" "${CYTO_MATH_LIB}" "${CYTO_BASE_LIB}"
    "${MEM3DG_LIB}"
    ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
    "${CURSES_LIBRARIES}"
    Threads::Threads
)

# so the executable finds libmem3dg at runtime without DYLD_/LD_LIBRARY_PATH
get_filename_component(_m3dg_libdir "${MEM3DG_LIB}" DIRECTORY)
set_target_properties(coupled_driver PROPERTIES
    BUILD_RPATH   "${_m3dg_libdir}"
    INSTALL_RPATH "${_m3dg_libdir}")

if(UNIX AND NOT APPLE)
  target_link_libraries(coupled_driver PRIVATE gfortran)
endif()
