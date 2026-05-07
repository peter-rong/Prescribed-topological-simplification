# Prescribed topological simplification

Reference code for **Persistence-guided Prescribed Topological Simplification** (SIGGRAPH 2026).

One executable supports **two usage modes** (same overall pipeline: filtration, then TopoMinCut min-cut step):

| Mode | When to use it |
|------|----------------|
| **Cubical** | Scalar field on a 3D **regular grid** (MRC volume) |
| **Tet** | Scalar field on an unstructured **tetrahedral mesh** (Gmsh `.msh` + per-vertex **alpha** file) |

TopoMinCut lives under [`TopoMinCut/`](TopoMinCut/) and is **linked in-process** (boundary matrix and alphas passed as Eigen structures; the main driver does not spawn an external TopoMinCut process).

## Dependencies

- **CMake** ≥ 3.10, **C++17**
- **Boost.Graph** — provide via your system, `CMAKE_PREFIX_PATH`, or [vcpkg](https://github.com/microsoft/vcpkg)
- **Eigen** — headers already included under [`TopoMinCut/dependencies`](TopoMinCut/dependencies)
- **PHAT** — vendored under [`TopoMinCut/phat`](TopoMinCut/phat)

**Platforms:** Developed and tested on **Windows** and **macOS**. **Linux** has not been exercised here; it should work in principle but may need small syntax or path tweaks.

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

## Usage

```text
PrescribedTopologicalSimplification <cubical|cubic|tet> ...arguments...
```

`-h`, `--help`, or `/?` prints a short reminder.

Run from a directory where relative **`output/`** is OK (both modes write artifacts there).

---

### Cubical mode

```text
PrescribedTopologicalSimplification cubical <input.mrc> <output_file> <output_file2> <adjustment> <dtype> [options]
```

Aliases: **`cubic`** (same as **`cubical`**).

| Argument / option | Description |
|-------------------|-------------|
| `input.mrc` | MRC scalar volume |
| `output_file`, `output_file2` | Basenames; files are placed under `output/` |
| `adjustment` | Offset added to the scalar field so the target isosurface is the **0-level set** of the volume |
| `dtype` | Dtype flag for the pipeline |
| `-a`, `--ascii` *(optional)* | ASCII IO mode |
| `--cpp_program <path>` *(optional)* | Ignored (backward-compatible scripts) |
| `--topK <n>` *(optional)* | TopoMinCut: skip longest `n` persistence pairs when building cuts |
| `--cavitySkip <n>` *(optional)* | TopoMinCut: skip the first `n` cavity (2D) features in generating sets |
| `--handleSkip <n>` *(optional)* | TopoMinCut: skip the first `n` handle (1D) features |
| `--componentSkip <n>` *(optional)* | TopoMinCut: skip the first `n` component (0D) features |

---

### Tet mode

```text
PrescribedTopologicalSimplification tet <mesh.msh> <alpha_file> <output_file> <output_file2> <adjustment> <dtype> [options]
```

| Argument / option | Description |
|-------------------|-------------|
| `mesh.msh` | Gmsh tet mesh |
| `alpha_file` | One scalar **alpha** per mesh vertex (e.g. signed distance to the surface); file format is consumed by `MSHReader::readVertexAlphas` |
| `output_file`, `output_file2` | Basenames under `output/` |
| `adjustment` | Same as cubical: offsets the field so the desired isosurface is the **0-level set** |
| `dtype` | Dtype flag for the pipeline |
| `-a`, `--ascii` *(optional)* | ASCII IO mode |
| `--cpp_program <path>` *(optional)* | Ignored (backward compatibility) |
| `--topK <n>` *(optional)* | Same as cubical |
| `--tet_labels <path>` *(optional)* | Tet label file |
| `--cavitySkip <n>` *(optional)* | Same as cubical (TopoMinCut) |
| `--handleSkip <n>` *(optional)* | Same as cubical |
| `--componentSkip <n>` *(optional)* | Same as cubical |

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
