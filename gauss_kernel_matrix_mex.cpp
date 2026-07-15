/*
 * gauss_kernel_matrix_mex.cpp
 * Precomputes the geometry-only half of the OPSCEA Gaussian proximity
 * kernel: the N_pts x N_elec matrix of exp(-distance^2/gsp) values between
 * every query point (brain vertex or slice pixel) and every electrode.
 *
 * Rationale: vertex/electrode positions are fixed for an entire OPSCEA
 * run - only the per-frame electrode weights change. gauss_kernel_mex
 * recomputed this distance/exp term from scratch on every single frame
 * even though it never changes. Computing it once here and combining it
 * with per-frame weights in plain MATLAB (a single matrix-vector multiply
 * for sum mode, or a multiply+max for max mode) eliminates that redundant
 * per-frame cost entirely.
 *
 * Usage:
 *   K = gauss_kernel_matrix_mex(pts, elecs, gsp)
 *
 * Inputs:
 *   pts   - N_pts x 3 double, query point coordinates (brain vertices or slice pixels)
 *   elecs - N_elec x 3 double, electrode coordinates
 *   gsp   - scalar double, Gaussian spread parameter
 *
 * Output:
 *   K - N_pts x N_elec double matrix, K(v,e) = exp(-dist(v,e)^2/gsp)
 *
 * Compile (macOS, no OpenMP):
 *   mex -O COPTIMFLAGS='-O3 -ffast-math' gauss_kernel_matrix_mex.cpp
 *
 * Compile with OpenMP (requires: brew install libomp):
 *   mex -O COPTIMFLAGS='-O3 -ffast-math -Xclang -fopenmp -I/usr/local/opt/libomp/include' ...
 *       LDFLAGS='$LDFLAGS -L/usr/local/opt/libomp/lib -lomp' gauss_kernel_matrix_mex.cpp
 */

#include "mex.h"
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 3)
        mexErrMsgIdAndTxt("gauss_kernel_matrix:nrhs",
            "3 inputs required: pts, elecs, gsp");
    if (nlhs > 1)
        mexErrMsgIdAndTxt("gauss_kernel_matrix:nlhs", "At most 1 output.");

    const double *pts   = mxGetPr(prhs[0]);
    const double *elecs = mxGetPr(prhs[1]);
    const double  gsp   = mxGetScalar(prhs[2]);

    const mwSize N_pts  = mxGetM(prhs[0]);
    const mwSize N_elec = mxGetM(prhs[1]);

    const double *px = pts;
    const double *py = pts + N_pts;
    const double *pz = pts + 2 * N_pts;

    const double *ex = elecs;
    const double *ey = elecs + N_elec;
    const double *ez = elecs + 2 * N_elec;

    const double neg_inv_gsp = -1.0 / gsp;

    plhs[0] = mxCreateDoubleMatrix(N_pts, N_elec, mxREAL);
    double *K = mxGetPr(plhs[0]);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (mwSize e = 0; e < N_elec; e++) {
        const double exx = ex[e], eyy = ey[e], ezz = ez[e];
        double *Kcol = K + e * N_pts; // column-major: column e starts at offset e*N_pts
        for (mwSize v = 0; v < N_pts; v++) {
            const double dx = px[v] - exx;
            const double dy = py[v] - eyy;
            const double dz = pz[v] - ezz;
            Kcol[v] = exp((dx*dx + dy*dy + dz*dz) * neg_inv_gsp);
        }
    }
}
