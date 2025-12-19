// laser_field.h
//
// This header defines functions for computing per‑vertex laser processing
// parameters (power, scan speed, pulse frequency) for the BFF mesh.
// The algorithm compensates for UV distortion, surface orientation, and
// curvature. It computes a scalar field on each mesh vertex and exports
// path samples as JSON. This module is designed to be integrated into
// Boundary First Flattening (BFF) without relying on OpenMesh.

#pragma once

#include <vector>
#include <string>
#include <memory>

// Forward declarations for BFF types. These are defined in the BFF
// library headers. We avoid including all of BFF here to keep this
// header lightweight. Consumers should include the appropriate BFF
// headers (e.g. <bff/mesh/Mesh.h>, <bff/project/Bff.h>) before
// including this file.

namespace bff {
class Bff;           // forward declaration of the main BFF project class
class Mesh;          // mesh type storing vertices, faces, corners
class Vertex;        // mesh vertex
class Corner;        // mesh corner (holds UV coordinates)
class Face;          // mesh face
class HalfEdge;      // halfedge structure
struct Vector;       // 2D or 3D vector (depending on context)
} // namespace bff

namespace laser {

// LaserParams encapsulates process settings used by the algorithm.
// Values are chosen based on typical laser marking systems. Users
// should tune these to match their hardware.
struct LaserParams {
    double targetAreaEnergy = 1.0;  // E_A*, J/mm^2: desired energy per unit area
    double spotSize        = 0.1;   // d, mm: nominal laser spot diameter
    double overlap         = 0.5;   // O, 0..1: pulse overlap ratio
    double Pmin            = 0.1;   // minimum power, W (1% rated)
    double Pmax            = 40.0;  // maximum power, W (40% rated)
    double fmin            = 10.0;  // minimum pulse frequency, Hz
    double fmax            = 200.0; // maximum pulse frequency, Hz
    double curvatureWeight = 0.0;   // γ: curvature compensation weight (optional)
    double beta            = 0.0;   // β: curvature influence on spot size
    double hatchUV         = 0.01;  // h_uv: UV hatch spacing (unitless)
    // Circular texture (dose) parameters in UV space
    bff::Vector circleCenter;       // (u0, v0) centre of high‑dose region
    double circleRadius    = 0.3;   // r: radius of high‑dose region in UV
    double doseInside      = 1.0;   // relative dose inside circle
    double doseOutside     = 0.5;   // relative dose outside circle

    // Path interpolation step.  When generating toolpath samples from
    // the vertex field, edges and faces will be sampled at roughly
    // this interval in millimetres.  A value of zero disables
    // interpolation and only vertex points are included.
    double interpolateStep = 0.0;
};

// VertexField stores computed laser parameters for a mesh vertex.
// index is the vertex index in the mesh; P, v, f store the computed
// power (W), scan speed (mm/s), and frequency (Hz), respectively.
struct VertexField {
    int    index;
    double P;
    double v;
    double f;
};

// Compute per‑vertex laser fields for the given BFF project and mesh.
//
// Arguments:
//  - bff: the active BFF project (already flattened). Its mesh holds
//    UV coordinates in each corner. The mesh is assumed to be triangulated.
//  - params: laser processing parameters (spot size, target energy, etc.).
//  - fieldsOut: output vector populated with one entry per vertex. The
//    order corresponds to vertex indices (fieldsOut[i] corresponds to
//    mesh->vertex[i]).
//
// The algorithm uses the following steps:
// 1. Compute per‑vertex normals and mean curvature using BFF geometry
//    utilities. Mean curvature is estimated via the discrete Laplace
//    operator (cotan weights). Normals and curvature are required for
//    angle and curvature compensation.
// 2. Compute UV distortion per vertex. For each face, compute the
//    ratio of its 3D area to its UV area; take the square root to
//    obtain an average length distortion factor. Accumulate these
//    distortions area‑weighted to vertices.
// 3. For each corner (or vertex), evaluate the texture dose field D
//    based on its UV coordinate: D = params.doseInside inside
//    the circle, or params.doseOutside outside. This modulates the
//    target energy.
// 4. For each vertex i, compute effective area energy E_A(i) =
//    params.targetAreaEnergy * D(i) * C_angle(i) * (1 + γ |H_i| d)
//    where C_angle compensates for surface tilt relative to the fixed
//    laser direction (Z axis) and γ |H| d adds optional curvature
//    compensation. If curvatureWeight is zero, curvature is ignored.
// 5. Convert the area energy E_A to a line energy E_L using the
//    surface hatch spacing h_s = params.hatchUV * s_L(i), where
//    s_L(i) is the length distortion factor from step 2. Then
//    E_L(i) = E_A(i) * h_s.
// 6. Determine a global pulse frequency f. The algorithm uses
//    params.fmax if possible; if the worst‑case required power at
//    fmax exceeds Pmax, scale down f accordingly but not below
//    params.fmin.
// 7. Compute base scan speed v_base = (1‑overlap) * spotSize * f.
//    For each vertex, compute preliminary power P_i = E_L(i) * v_base.
//    Clamp P_i to [Pmin, Pmax]; adjust v_i accordingly if P was
//    clamped (slower if underpowered). Set f_i = f for all vertices.
//
// Returns true on success; false if an error occurs (e.g. missing UV).
bool computeLaserFields(bff::Bff* project,
                        const LaserParams& params,
                        std::vector<VertexField>& fieldsOut);

// Export the computed vertex field to a JSON file. This helper
// writes a JSON array with one object per vertex containing the
// vertex index, coordinates, and computed parameters. The mesh must
// still be accessible via the BFF project to query positions.
//
// Example JSON structure:
// [
//   {"index":0,"x":...,"y":...,"z":...,"P":...,"v":...,"f":...},
//   ...
// ]
bool exportLaserFieldAsJson(bff::Bff* project,
                            const std::vector<VertexField>& fields,
                            const std::string& filepath);

// A single point on the toolpath.  Each PathPoint records the 3D
// coordinates (x,y,z) along with the interpolated laser parameters
// (P, v, f).  Path points are generated by subdividing faces and
// edges of the mesh to create a continuous sampling of the vertex
// field.  A separate path point can be used by downstream code to
// construct scan lines or motion commands.
struct PathPoint {
    double x, y, z;
    double P;
    double v;
    double f;
};

// Interpolate the vertex field onto a set of path points.  This
// routine samples each triangle in the mesh at a resolution
// determined by step (approximate spacing in mm).  For each
// subdivision point, the barycentric coordinates within the face are
// used to blend the three vertex parameters.  The resulting points
// are appended to the output vector.  If step <= 0, no sampling is
// performed and only the original vertex positions are converted to
// path points.
void interpolateFieldToPath(bff::Bff* project,
                            const std::vector<VertexField>& fields,
                            double step,
                            std::vector<PathPoint>& outPath);

// Export a path (generated by interpolateFieldToPath) as a JSON
// array.  Each entry contains x, y, z, P, v, and f.  Returns true
// on success.
bool exportPathAsJson(const std::vector<PathPoint>& path,
                      const std::string& filepath);

} // namespace laser