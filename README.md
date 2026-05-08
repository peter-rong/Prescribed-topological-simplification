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
PrescribedTopologicalSimplification cubical <input.mrc> <adjustment> [options]
```

Aliases: **`cubic`** (same as **`cubical`**).

| Argument / option | Description |
|-------------------|-------------|
| `input.mrc` | MRC scalar volume |
| `adjustment` | Offset added to the scalar field so the target isosurface is the **0-level set** of the volume |
| `-a`, `--ascii` *(optional)* | ASCII IO mode |
| `--topK <n>` *(optional)* | Ignore the *K* most persistent topological features. You do not need this flag if you set any of `--cavitySkip`, `--handleSkip`, or `--componentSkip`. |
| `--cavitySkip <n>` *(optional)* | How many cavity (2D) features to preserve |
| `--handleSkip <n>` *(optional)* | How many handle (1D) features to preserve |
| `--componentSkip <n>` *(optional)* | How many component (0D) features to preserve |

Cubical mode also writes a filtered complex dump and per-simplex alphas to fixed paths: **`output/cubical_filtered_simplices.txt`** and **`output/cubical_filtered_alphas.txt`**.

---

### Tet mode

```text
PrescribedTopologicalSimplification tet <mesh.msh> <alpha_file> <adjustment> [options]
```

| Argument / option | Description |
|-------------------|-------------|
| `mesh.msh` | Gmsh tet mesh |
| `alpha_file` | One scalar **alpha** per mesh vertex (e.g. signed distance to the surface); file format is consumed by `MSHReader::readVertexAlphas` |
| `adjustment` | Offset added to the scalar field so the target isosurface is the **0-level set** of the volume |
| `-a`, `--ascii` *(optional)* | ASCII IO mode |
| `--topK <n>` *(optional)* | Ignore the *K* most persistent topological features. You do not need this flag if you set any of `--cavitySkip`, `--handleSkip`, or `--componentSkip`. |
| `--tet_labels <path>` *(optional)* | Optional Tet label file to perform subdivision on input mesh for a more accurate input representation|
| `--cavitySkip <n>` *(optional)* | How many cavity (2D) features to preserve |
| `--handleSkip <n>` *(optional)* | How many handle (1D) features to preserve |
| `--componentSkip <n>` *(optional)* | How many component (0D) features to preserve |

---

## Layout

| Path | Role |
|------|------|
| `TopoMinCut\` | Core code for topological simplification |
| `main.cpp` | Mode dispatch |
| `cubical_mode.*` | Grid / MRC |
| `tet_mode.*` | Tet mesh / `.msh` |
| `marching_cubes.*` | Isosurface (cubical path) |
| `mrc_reader.h`, `msh_reader.h` | Readers |

## Citation

If you use this code academically, cite the paper *(add bibliographic entry when available)*.
