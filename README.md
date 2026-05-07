# Prescribed topological simplification

Reference code for **Persistence-guided Prescribed Topological Simplification** (SIGGRAPH 2026).

One executable supports **two usage modes** (same overall pipeline: filtration, then TopoMinCut min-cut step):

| Mode | When to use it |
|------|----------------|
| **Cubical** | Scalar field on a 3D **regular grid** (MRC volume) |
| **Tet** | Scalar field on an unstructured **tetrahedral mesh** (Gmsh `.msh` + boundary file) |

TopoMinCut lives under [`TopoMinCut/`](TopoMinCut/) and is **linked in-process** (boundary matrix and alphas passed as Eigen structures; the main driver does not spawn an external TopoMinCut process).

## Dependencies

- **CMake** ≥ 3.10, **C++17**
- **Boost.Graph** — provide via your system, `CMAKE_PREFIX_PATH`, or [vcpkg](https://github.com/microsoft/vcpkg)
- **Eigen** — headers already included under [`TopoMinCut/dependencies`](TopoMinCut/dependencies)
- **PHAT** — vendored under [`TopoMinCut/phat`](TopoMinCut/phat)

## Build

**Windows** (Boost via vcpkg):

```bash
cmake -B build -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build --config Release
```

**Other** (if Boost is already discoverable by CMake):

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Output binary:

- Linux/macOS: `build/PrescribedTopologicalSimplification`
- Windows (VS generator): `build/Release/PrescribedTopologicalSimplification.exe`

### Tet mode: parallel `std::sort` (optional)

If linking fails on some Linux setups (e.g. missing TBB for `<execution>`), uncomment in root `CMakeLists.txt`:

```cmake
target_compile_definitions(PrescribedTopologicalSimplification PRIVATE TOPOTET_NO_PAR_SORT)
```

## Command line

```text
PrescribedTopologicalSimplification <cubical|cubic|tet> ...arguments...
```

`-h`, `--help`, or `/?` prints a short reminder.

Run from a directory where relative **`output/`** is OK (both modes write artifacts there).

---

### Cubical mode

```text
PrescribedTopologicalSimplification cubical <input.mrc> <output_file> <output_file2> <adjustment> <dtype> [optional flags...]
```

Aliases: **`cubic`** (same as **`cubical`**).

**Required**

| Argument | Description |
|----------|-------------|
| `input.mrc` | MRC scalar volume |
| `output_file`, `output_file2` | Basenames; files are placed under `output/` |
| `adjustment` | Scalar used in the pipeline (with `core` / `neighborhood`) |
| `dtype` | Dtype flag for the pipeline |

**Optional flags**

| Flag | Description |
|------|-------------|
| `-a`, `--ascii` | ASCII IO mode |
| `--cpp_program <path>` | Ignored (kept for backward-compatible scripts) |
| `--topK <n>` | TopoMinCut: skip longest `n` persistence pairs when building cuts |
| `--core <value>` | Core threshold (see implementation after adjustment) |
| `--neighborhood <value>` | Neighborhood threshold (after adjustment) |

---

### Tet mode

```text
PrescribedTopologicalSimplification tet <mesh.msh> <boundary_file> <output_file> <output_file2> <adjustment> <dtype> [optional flags...]
```

**Required**

| Argument | Description |
|----------|-------------|
| `mesh.msh` | Gmsh tet mesh |
| `boundary_file` | Boundary definition |
| `output_file`, `output_file2` | Basenames under `output/` |
| `adjustment`, `dtype` | Same roles as cubical mode |

**Optional flags**

| Flag | Description |
|------|-------------|
| `-a`, `--ascii` | ASCII IO mode |
| `--cpp_program <path>` | Ignored (backward compatibility) |
| `--topK <n>` | Same as cubical |
| `--core <value>` | Same as cubical |
| `--neighborhood <value>` | Same as cubical |
| `--tet_labels <path>` | Optional tet label file |
| `--tetMetricsLog <jsonl>` | Append driver timing (JSON Lines) |
| `--cavitySkip <n>` | Skip first `n` cavity features (per driver logic) |
| `--handleSkip <n>` | Skip first `n` handle features |
| `--componentSkip <n>` | Skip first `n` component features |
| `--topoMinCutMetricsLog <path>` | Accepted for compatibility; TopoMinCut no longer writes this log |

---

### TopoMinCut as a library

For custom integrations, see [`TopoMinCut/include/topomincut_runner.hpp`](TopoMinCut/include/topomincut_runner.hpp) (`runFromEigenSparse`, `runFromBoundaryMatrix`). A small standalone **`TopoMinCut`** executable is still built for file-based matrix/alpha experiments; the main prescribed simplification driver does not require it.

---

## Layout

| Path | Role |
|------|------|
| `main.cpp` | Mode dispatch |
| `cubical_mode.*` | Grid / MRC |
| `tet_mode.*` | Tet mesh / `.msh` |
| `marching_cubes.*` | Isosurface (cubical path) |
| `mrc_reader.h`, `msh_reader.h` | Readers |

## Citation

If you use this code academically, cite the paper *(add bibliographic entry when available)*.
