#pragma once

namespace prescribed_topo {
namespace tet {

/** Tet mesh mode: Gmsh .msh + per-vertex alpha file (same CLI shape as legacy TopoTet after the mode token). */
int runTetMode(int argc, char* argv[]);

} // namespace tet
} // namespace prescribed_topo
