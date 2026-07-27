# Building the coupled Mem3DG + Cytosim driver

Three separate builds. They stay separate on purpose: Mem3DG is C++20 and
exposes Eigen and geometry-central in its headers, Cytosim is C++11. The only
thing that crosses between them is `<mem3dg/capi/mem3dg_c.h>`, a plain-C ABI
over an opaque handle, so neither project has to compile the other's headers.

Set these once:

```bash
export MEM3DG=/Users/hakenuwa/Desktop/Mem3DG
export CYTOSIM=/Users/hakenuwa/Desktop/cytosim-lee-lab-2
export BUILD=/Users/hakenuwa/Desktop/Mem3DG/coupled-build          # anywhere you like
```

## 1. libmem3dg (C++20)

```bash
cmake -S "$MEM3DG" -B "$BUILD/mem3dg" -DCMAKE_BUILD_TYPE=Release -DBUILD_PYMEM3DG=OFF -DBUILD_C_API=ON && cmake --build "$BUILD/mem3dg" -j8
```

`BUILD_C_API` is ON by default; it compiles `src/capi/mem3dg_c.cpp` into the
shared library, which is the single place C linkage meets C++20. `BUILD_PYMEM3DG=OFF`
just skips the Python extension, which the driver does not use.

Produces `$BUILD/mem3dg/lib/libmem3dg.dylib` (`.so` on Linux). Check the ABI is
actually exported:

```bash
nm -gU "$BUILD/mem3dg/lib/libmem3dg.dylib" | grep -c mem3dg_system_
```

Expect 27. If it is 0, the symbols were hidden — the repo default visibility is
hidden, so every entry point needs `DLL_PUBLIC`.

## 2. Cytosim static libraries (C++11)

```bash
cmake -S "$CYTOSIM" -B "$BUILD/cytosim" -DCMAKE_BUILD_TYPE=Release -DMEM3DG_ROOT="$MEM3DG" -DMAKE_SIM=ON -DMAKE_PLAY=OFF -DMAKE_TOOLS=OFF -DMAKE_TESTS=OFF && cmake --build "$BUILD/cytosim" --target cytosim cytomath cytobase -j8
```

**`-DMEM3DG_ROOT` is what enables the `mem3dg` space shape.** Without it,
`src/sim/CMakeLists.txt` skips the include directory and `HAS_MEM3DG`, and any
config using `shape = mem3dg` dies at parse time with

```
unknown space:shape `mem3dg'
```

Only the *header* is needed here, not the library: Cytosim is a static archive,
so the `mem3dg_*` symbols stay unresolved until the driver links. That means
this step does not depend on step 1 having finished — only on the header
existing in the source tree.

Produces `$BUILD/cytosim/tar/libcytosim.a`, `libcytomath.a`, `libcytobase.a`.
Verify the space really compiled in:

```bash
ar t "$BUILD/cytosim/tar/libcytosim.a" | grep mem3dg
```

Expect `space_mem3dg.cc.o`. If it is missing, the registration in
`space_prop.cc` / `space_set.cc` is not being compiled and you will get the
`unknown space:shape` error again.

## 3. The driver (C++11)

```bash
cmake -S "$MEM3DG/coupled_driver" -B "$BUILD/driver" -DCMAKE_BUILD_TYPE=Release -DCYTOSIM_DIR="$CYTOSIM" -DCYTOSIM_BUILD="$BUILD/cytosim" -DMEM3DG_ROOT="$MEM3DG" -DMEM3DG_LIB="$BUILD/mem3dg/lib/libmem3dg.dylib" && cmake --build "$BUILD/driver" -j8
```

`CYTOSIM_DIR` is the source root (headers), `CYTOSIM_BUILD` is the build dir
(the `tar/` archives) — they are different paths and both are required. The
CMakeLists sets `BUILD_RPATH` to the libmem3dg directory, so the executable
finds it at runtime without `DYLD_LIBRARY_PATH`.

Produces `$BUILD/driver/coupled_driver`.

## DIM must match

Cytosim bakes the spatial dimension in at compile time. The driver defaults to
`-DDIMENSION=3` to match Mem3DG's 3D meshes; if you build Cytosim for another
dimension, pass the same value to the driver or the objects will disagree about
`Vector`'s size and you get silent memory corruption, not a link error.

## Smoke test

```bash
cd "$MEM3DG/tests/test-membrane-surface" && "$BUILD/driver/coupled_driver" config-nucleate.cym ../icosphere.obj 300 0.0035 10 1 100 /tmp/smoke 1 8.22e-6 100.0 100.0 0.98 0.0 0
```

Healthy output has `|F_fib|` equal to `|F_mem|` on every line (action-reaction),
and writes `membrane_%04d.obj`, `fibers_%04d.txt`, `points_%04d.txt` into
`/tmp/smoke`.

## Driver arguments

```
coupled_driver <config> <mesh> <nSteps> <shell/cutoff> <stiffness> <strang>
               <saveEvery> <outDir> <mode> <Kb> <Ksg> <Kv> <areaFrac> <Kse>
               <remesh> <maxEdge> <remeshEvery>
```

`mode` 0 = symmetric repulsion (works on open meshes, does not confine),
1 = one-sided wall (confines filaments inside, needs a closed mesh).
`strang` 1 = Strang splitting (error ~dt²), 0 = Lie–Trotter (~dt).
