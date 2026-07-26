# ============================================================================
# UseMem3DGInCytosim.cmake
#
# Injected into Cytosim's configure step WITHOUT editing any Cytosim file:
#
#   cmake -S /path/to/cytosim-lee-lab-2 -B build-cyto \
#     -DCMAKE_PROJECT_cytosim_INCLUDE=/Users/hakenuwa/Desktop/Mem3DG/cmake/UseMem3DGInCytosim.cmake \
#     -DMEM3DG_ROOT=/Users/hakenuwa/Desktop/Mem3DG \
#     -DMEM3DG_LIB=/Users/hakenuwa/Desktop/Mem3DG/build/lib/libmem3dg.dylib
#
# CMake includes this file immediately after Cytosim's `project("cytosim" ...)`
# call — which is BEFORE `add_subdirectory("src")` where sim/play/tools are
# defined.  The link_libraries()/include_directories() calls below therefore
# apply globally to every Cytosim target created afterward, so `sim`, `play`,
# etc. all link libmem3dg and can #include <mem3dg/capi/mem3dg_c.h>.
# ============================================================================

if(NOT DEFINED MEM3DG_ROOT)
  message(FATAL_ERROR "Set -DMEM3DG_ROOT=/abs/path/to/Mem3DG")
endif()
if(NOT DEFINED MEM3DG_LIB)
  # Default guess: the in-tree build directory.  CMAKE_SHARED_LIBRARY_SUFFIX is
  # ".dylib" on macOS and ".so" on Linux, so one expression covers both.
  set(MEM3DG_LIB
      "${MEM3DG_ROOT}/build/lib/libmem3dg${CMAKE_SHARED_LIBRARY_SUFFIX}")
endif()
if(NOT EXISTS "${MEM3DG_LIB}")
  message(FATAL_ERROR "libmem3dg not found at ${MEM3DG_LIB}. "
                      "Build it first:  cmake --build ${MEM3DG_ROOT}/build --target mem3dg")
endif()

message(STATUS "[mem3dg] linking Cytosim targets against ${MEM3DG_LIB}")

add_library(mem3dg_shared SHARED IMPORTED GLOBAL)
set_target_properties(mem3dg_shared PROPERTIES
  IMPORTED_LOCATION            "${MEM3DG_LIB}"
  INTERFACE_INCLUDE_DIRECTORIES "${MEM3DG_ROOT}/include")

# Applies to every target defined after this point (i.e. all of Cytosim's).
link_libraries(mem3dg_shared)

# So <mem3dg/capi/mem3dg_c.h> resolves in any Cytosim .cc — no per-target edit.
include_directories("${MEM3DG_ROOT}/include")

# Embed the library's directory in the executables' RPATH so they find
# libmem3dg at runtime without needing LD_LIBRARY_PATH / DYLD_LIBRARY_PATH.
get_filename_component(_mem3dg_libdir "${MEM3DG_LIB}" DIRECTORY)
list(APPEND CMAKE_BUILD_RPATH "${_mem3dg_libdir}")
list(APPEND CMAKE_INSTALL_RPATH "${_mem3dg_libdir}")
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
