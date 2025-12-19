// laser_field.cpp
//
// Implementation of the per‑vertex laser processing parameter computation
// for the Boundary First Flattening project. See laser_field.h for a
// description of the algorithm.

#include "laser_field.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>

// Include BFF core headers. Consumers should ensure these headers
// are available in the include path when compiling this file with
// BFF. We include Mesh and geometry utilities to access mesh
// topology and compute geometric quantities.
#include "bff/mesh/Mesh.h"
#include "bff/mesh/Vertex.h"
#include "bff/mesh/Corner.h"
#include "bff/mesh/Face.h"
#include "bff/mesh/HalfEdge.h"
#include "bff/mesh/GeometryUtils.h"
#include "bff/project/Bff.h"

namespace laser {

// Helper function to compute length distortion factor s_L per vertex.
static std::vector<double> computeLengthDistortion(bff::Mesh& mesh) {
    size_t nV = mesh.vertices.size();
    std::vector<double> sL(nV, 1.0);
    // For each vertex, accumulate weighted s_L contributions
    // Use the method: s_L(i) = sum_f sqrt(area3D/areaUV) * area3D / sum_f area3D
    // If areaUV is zero or negative, ignore that face.
    for (auto& v : mesh.vertices) {
        double areaSum = 0.0;
        double slSum   = 0.0;
        // iterate over faces adjacent to v using halfedge fan
        auto hStart = v.halfEdge();
        auto h = hStart;
        if (!h) continue;
        do {
            auto f = h->face();
            if (f && f->isReal()) {
                double area3D = bff::normal(f, false).norm() * 0.5;
                double areaUV = bff::areaUV(f);
                if (areaUV > 1e-12) {
                    double sLface = std::sqrt(area3D / areaUV);
                    areaSum += area3D;
                    slSum   += sLface * area3D;
                }
            }
            // move to next face around the vertex
            if (!h->onBoundary) {
                h = h->flip()->next();
            } else {
                // boundary: still go to the next but stop after one loop
                h = h->flip();
                if (!h->onBoundary) h = h->next();
                else break;
            }
        } while (h != hStart && h);
        if (areaSum > 1e-12) sL[v.index] = slSum / areaSum;
    }
    return sL;
}

// Helper function to compute mean curvature magnitude per vertex using
// the discrete Laplace operator (cotan weights). For simplicity,
// curvatureWeight defaults to zero and can be ignored. If curvatureWeight
// is non‑zero, this function returns |H| values, else zeros.
static std::vector<double> computeMeanCurvature(bff::Mesh& mesh) {
    size_t nV = mesh.vertices.size();
    std::vector<double> H(nV, 0.0);
    // Precompute per vertex area (one third of adjacent face areas)
    std::vector<double> area(nV, 0.0);
    for (auto& f : mesh.faces) {
        if (!f.isReal()) continue;
        double a = bff::area(&f);
        // accumulate to each corner's vertex
        auto h = f.halfEdge();
        for (int i = 0; i < 3; ++i) {
            auto v = h->vertex();
            area[v->index] += a / 3.0;
            h = h->next();
        }
    }
    // Compute Laplacian vector and curvature magnitude
    for (auto& v : mesh.vertices) {
        int i = v.index;
        bff::Vector laplace(0.0, 0.0, 0.0);
        double Ai = area[i];
        if (Ai < 1e-12) continue;
        auto hStart = v.halfEdge();
        auto h = hStart;
        if (!h) continue;
        do {
            // compute cotan weights for edge (v, vj) and (vj,v)
            auto hj = h;
            // weight for halfedge hj
            double cotanSum = bff::cotan(hj);
            // add weight of opposite face (flip)
            auto hf = hj->flip();
            if (hf && !hf->onBoundary) {
                cotanSum += bff::cotan(hf);
            }
            // vector difference to neighbor vertex
            bff::Vector diff = hj->next()->vertex()->position - v.position;
            laplace += cotanSum * diff;
            // next halfedge around vertex
            if (!hj->onBoundary) {
                h = hj->flip()->next();
            } else {
                h = hj->flip();
                if (!h->onBoundary) h = h->next();
                else break;
            }
        } while (h != hStart && h);
        // mean curvature normal = 1/(2A) * laplace
        laplace *= (1.0 / (2.0 * Ai));
        double Hmag = 0.5 * std::sqrt(laplace.dot(laplace));
        H[i] = Hmag;
    }
    return H;
}

bool computeLaserFields(bff::Bff* project,
                        const LaserParams& params,
                        std::vector<VertexField>& fieldsOut) {
    if (!project) return false;
    // Access underlying mesh
    bff::Mesh& mesh = project->mesh;
    size_t nV = mesh.vertices.size();
    fieldsOut.clear();
    fieldsOut.resize(nV);
    // Precompute per‑vertex normal (for angle compensation)
    std::vector<bff::Vector> normals(nV);
    for (auto& v : mesh.vertices) {
        normals[v.index] = bff::normal(&v);
    }
    // Compute length distortion factor s_L per vertex
    std::vector<double> sL = computeLengthDistortion(mesh);
    // Compute curvature magnitude per vertex (optional)
    std::vector<double> H = params.curvatureWeight != 0.0 ? computeMeanCurvature(mesh) : std::vector<double>(nV, 0.0);
    // Determine global pulse frequency. Start with fmax and compute
    // worst required power at baseline speed.
    double f = params.fmax;
    // We'll compute worst required P at fmax and adjust f if needed.
    auto computeWorstPower = [&](double freq) -> double {
        double v_base = (1.0 - params.overlap) * params.spotSize * freq;
        double worstP = 0.0;
        for (auto& v : mesh.vertices) {
            int i = v.index;
            // get first wedge (corner) to sample uv
            auto w = v.wedge();
            if (!w) continue;
            bff::Vector uv = w->uv;
            // compute dose
            double dx = uv.x - params.circleCenter.x;
            double dy = uv.y - params.circleCenter.y;
            double dist2 = dx*dx + dy*dy;
            double D = (dist2 <= params.circleRadius * params.circleRadius) ? params.doseInside : params.doseOutside;
            // angle compensation
            double cosT = std::max(1e-6, normals[i].z);
            double C_angle = 1.0 / cosT;
            double C_curv  = 1.0 + params.curvatureWeight * H[i] * params.spotSize;
            double E_A     = params.targetAreaEnergy * D * C_angle * C_curv;
            double h_s     = params.hatchUV * sL[i];
            double E_L     = E_A * h_s;
            double P_req   = E_L * v_base;
            if (P_req > worstP) worstP = P_req;
        }
        return worstP;
    };
    // Evaluate worst power at fmax
    double worstPower = computeWorstPower(f);
    if (worstPower > params.Pmax) {
        // reduce f proportionally to fit Pmax
        double scale = params.Pmax / worstPower;
        double fNew = f * scale;
        if (fNew < params.fmin) fNew = params.fmin;
        f = fNew;
    }
    // Compute base speed with final f
    double v_base = (1.0 - params.overlap) * params.spotSize * f;
    // Loop through vertices to compute final fields
    for (auto& v : mesh.vertices) {
        int i = v.index;
        VertexField vf; vf.index = i; vf.f = f;
        auto w = v.wedge();
        if (!w) {
            // if no wedge, set defaults
            vf.P = params.Pmin;
            vf.v = v_base;
            fieldsOut[i] = vf;
            continue;
        }
        bff::Vector uv = w->uv;
        double dx = uv.x - params.circleCenter.x;
        double dy = uv.y - params.circleCenter.y;
        double dist2 = dx*dx + dy*dy;
        double D = (dist2 <= params.circleRadius * params.circleRadius) ? params.doseInside : params.doseOutside;
        double cosT = std::max(1e-6, normals[i].z);
        double C_angle = 1.0 / cosT;
        double C_curv  = 1.0 + params.curvatureWeight * H[i] * params.spotSize;
        double E_A     = params.targetAreaEnergy * D * C_angle * C_curv;
        double h_s     = params.hatchUV * sL[i];
        double E_L     = E_A * h_s;
        double P_req   = E_L * v_base;
        // clamp power
        double P = P_req;
        if (P > params.Pmax) P = params.Pmax;
        if (P < params.Pmin) P = params.Pmin;
        // adjust speed if underpowered (P was clamped down)
        double v_i = v_base;
        if (P < P_req) {
            // slow down to satisfy energy density
            v_i = P / (E_L + 1e-12);
        }
        // we choose not to speed up if P > P_req (overpowering is acceptable)
        vf.P = P;
        vf.v = v_i;
        fieldsOut[i] = vf;
    }
    return true;
}

bool exportLaserFieldAsJson(bff::Bff* project,
                            const std::vector<VertexField>& fields,
                            const std::string& filepath) {
    if (!project) return false;
    bff::Mesh& mesh = project->mesh;
    std::ostringstream ss;
    ss << "[\n";
    bool first = true;
    for (const auto& vf : fields) {
        int idx = vf.index;
        if (idx < 0 || static_cast<size_t>(idx) >= mesh.vertices.size()) {
            continue;
        }
        const bff::Vertex& v = mesh.vertices[idx];
        if (!first) ss << ",\n";
        ss << "  {\"index\":" << idx
            << ",\"x\":" << v.position.x
            << ",\"y\":" << v.position.y
            << ",\"z\":" << v.position.z
            << ",\"P\":" << vf.P
            << ",\"v\":" << vf.v
            << ",\"f\":" << vf.f
            << "}";
        first = false;
    }
    std::ofstream ofs(filepath);
    if (!ofs) return false;
    ofs << ss.str() << "\n]";
    ofs.close();
    return true;
}

void interpolateFieldToPath(bff::Bff* project,
                            const std::vector<VertexField>& fields,
                            double step,
                            std::vector<PathPoint>& outPath) {
    outPath.clear();
    if (!project) return;
    bff::Mesh& mesh = project->mesh;
    size_t nV = mesh.vertices.size();
    // create a lookup by vertex index to field values
    std::vector<const VertexField*> vfield(nV, nullptr);
    for (const auto& vf : fields) {
        if (vf.index >= 0 && static_cast<size_t>(vf.index) < nV) {
            vfield[vf.index] = &vf;
        }
    }
    // If no step or non‑positive, simply output each vertex
    if (step <= 0.0) {
        for (const auto& v : mesh.vertices) {
            int i = v.index;
            const VertexField* vf = vfield[i];
            if (!vf) continue;
            PathPoint pp;
            pp.x = v.position.x;
            pp.y = v.position.y;
            pp.z = v.position.z;
            pp.P = vf->P;
            pp.v = vf->v;
            pp.f = vf->f;
            outPath.push_back(pp);
        }
        return;
    }
    // sample each face by subdividing into approximately equal segments along edges
    for (const auto& f : mesh.faces) {
        if (!f.isReal()) continue;
        // get vertices and fields
        std::array<bff::Vertex*,3> verts;
        std::array<const VertexField*,3> fieldsV;
        auto he = f.halfEdge();
        for (int k = 0; k < 3; ++k) {
            bff::Vertex* v = he->vertex();
            verts[k] = v;
            fieldsV[k] = vfield[v->index];
            he = he->next();
        }
        // if any field missing skip face
        if (!fieldsV[0] || !fieldsV[1] || !fieldsV[2]) continue;
        // compute edge lengths using vector norms
        bff::Vector d01 = verts[1]->position - verts[0]->position;
        bff::Vector d12 = verts[2]->position - verts[1]->position;
        bff::Vector d20 = verts[0]->position - verts[2]->position;
        double len0 = d01.norm();
        double len1 = d12.norm();
        double len2 = d20.norm();
        double maxLen = std::max(len0, std::max(len1, len2));
        int nSteps = static_cast<int>(std::ceil(maxLen / step));
        if (nSteps < 1) nSteps = 1;
        // sample barycentric coordinates u,v such that u,v >=0, u+v<=1
        // we create a simple grid over the triangle
        for (int iu = 0; iu <= nSteps; ++iu) {
            for (int iv = 0; iv <= nSteps - iu; ++iv) {
                double a = static_cast<double>(iu) / nSteps;
                double b = static_cast<double>(iv) / nSteps;
                double c = 1.0 - a - b;
                // barycentric combination of vertex positions
                double px = c*verts[0]->position.x + a*verts[1]->position.x + b*verts[2]->position.x;
                double py = c*verts[0]->position.y + a*verts[1]->position.y + b*verts[2]->position.y;
                double pz = c*verts[0]->position.z + a*verts[1]->position.z + b*verts[2]->position.z;
                // interpolate fields linearly
                double P  = c*fieldsV[0]->P + a*fieldsV[1]->P + b*fieldsV[2]->P;
                double v_i = c*fieldsV[0]->v + a*fieldsV[1]->v + b*fieldsV[2]->v;
                double f_i = c*fieldsV[0]->f + a*fieldsV[1]->f + b*fieldsV[2]->f;
                PathPoint pp;
                pp.x = px; pp.y = py; pp.z = pz;
                pp.P = P; pp.v = v_i; pp.f = f_i;
                outPath.push_back(pp);
            }
        }
    }
}

bool exportPathAsJson(const std::vector<PathPoint>& path,
                      const std::string& filepath) {
    std::ofstream ofs(filepath);
    if (!ofs) return false;

    ofs << "[\n";
    for (size_t i = 0; i < path.size(); ++i) {
        const auto& p = path[i];
        ofs << "  {\"x\":" << p.x
            << ",\"y\":" << p.y
            << ",\"z\":" << p.z
            << ",\"P\":" << p.P
            << ",\"v\":" << p.v
            << ",\"f\":" << p.f
            << "}";
        if (i + 1 < path.size()) ofs << ",";
        ofs << "\n";
    }
    ofs << "]";
    return true;
}
} // namespace laser
