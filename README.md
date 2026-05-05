# Prescribed topological simplification

Reference implementation for **Persistence-based Prescribed Topological Simplification** in two geometry backends:

| Mode | Domain | Primary input |
|------|--------|----------------|
| **Cubical** | Structured 3D grid | MRC volume (scalar field on voxels) |
| **Tet** | Unstructured tetrahedral mesh | Gmsh `.msh` plus a boundary-definition file |

Both modes share the same high-level pipeline (cell complex, filtration, and calls into an external min-cut solver) but use **different input parsers and combinatorial complexes**, chosen explicitly at the command line.

## Requirements

- **CMake** ≥ 3.10  
- **C++17** compiler (GCC, Clang, or MSVC)

Optional:

- **TopoMinCut** (or compatible driver): both modes invoke an external program passed via `--cpp_program`. Set this to your built binary so runs are reproducible across machines.

## Build

Configure and build from this directory:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

The executable is produced as:

- **Linux/macOS:** `build/PrescribedTopologicalSimplification`  
- **Windows (multi-config generators):** `build/Release/PrescribedTopologicalSimplification.exe`

### Parallel sorting (tet mode, optional)

Tet mode may use `std::sort` with `<execution>` when available. If linking fails (for example, missing TBB on some Linux setups), uncomment in `CMakeLists.txt`:

```cmake
target_compile_definitions(PrescribedTopologicalSimplification PRIVATE TOPOTET_NO_PAR_SORT)
```

## Usage overview

The first argument selects the backend; remaining arguments match the legacy **TopoCubical** / **TopoTet** CLIs (without repeating the executable name).

```text
PrescribedTopologicalSimplification <cubical|tet> ...
```

Use **`--help`** / **`-h`** / **`/?`** for a short reminder printed by the binary.

**Working directory:** many artifacts are written under a relative `output/` folder (meshes, skeleton PLY files, logs). Run the tool from the project directory where you want that `output/` tree created, or ensure relative paths resolve as intended.

---

## Cubical mode (`cubical` or `cubic`)

For **regular-grid** data read from an **MRC** file.

```text
PrescribedTopologicalSimplification cubical <input.mrc> <output_file> <output_file2> <adjustment> <dtype> [options]
```

| Argument | Meaning |
|----------|---------|
| `input.mrc` | Path to the scalar volume (MRC format). |
| `output_file`, `output_file2` | Basenames for intermediate matrices / alpha files (prepended with `output/` internally). |
| `adjustment` | Scalar shift applied in pipeline (see implementation / paper). |
| `dtype` | Dtype flag consumed by the pipeline (same as original cubical driver). |

**Optional flags**

| Flag | Purpose |
|------|---------|
| `-a`, `--ascii` | ASCII-related IO mode (same semantics as original). |
| `--cpp_program <path>` | Path to **TopoMinCut** (or compatible) executable. |
| `--topK <n>` | Passed through to the external solver. |
| `--core <value>` | Core threshold (after adjustment in driver). |
| `--neighborhood <value>` | Neighborhood threshold (after adjustment in driver). |

---

## Tet mode (`tet`)

For **tetrahedral meshes** from **Gmsh** plus boundary information.

```text
PrescribedTopologicalSimplification tet <mesh.msh> <boundary_file> <output_file> <output_file2> <adjustment> <dtype> [options]
```

| Argument | Meaning |
|----------|---------|
| `mesh.msh` | Tet mesh (Gmsh format). |
| `boundary_file` | Boundary data expected by the original tet driver. |
| `output_file`, `output_file2` | Same role as cubical mode (under `output/`). |
| `adjustment`, `dtype` | Same convention as the original tet driver. |

**Optional flags**

Includes the same common flags as cubical mode (`-a`, `--cpp_program`, `--topK`, `--core`, `--neighborhood`), plus:

| Flag | Purpose |
|------|---------|
| `--tet_labels <path>` | Optional tet label file. |
| `--tetMetricsLog <jsonl>` | Append timing metrics (JSON Lines). |
| `--cavitySkip <n>` | Skip cavities (driver-specific). |
| `--handleSkip <n>` | Skip handles. |
| `--componentSkip <n>` | Skip components. |
| `--topoMinCutMetricsLog <jsonl>` | Metrics log for TopoMinCut. |

---

## Source layout

| File | Role |
|------|------|
| `main.cpp` | Mode dispatch and CLI parsing. |
| `cubical_mode.cpp` / `cubical_mode.h` | Grid / cubical complex and MRC IO. |
| `tet_mode.cpp` / `tet_mode.h` | Tet complex and `.msh` IO. |
| `marching_cubes.cpp`, `marching_cubes.h` | Isosurface extraction (cubical path). |
| `mrc_reader.h` | MRC reader (header-only implementation). |
| `msh_reader.h` | Gmsh reader (header-only implementation). |

---

## Citation

If you use this code in academic work, please cite the accompanying paper *(add citation here once available)*.
