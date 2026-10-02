#ifndef LAPLACE_H_
#define LAPLACE_H_

#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

#include <stdint.h>

static inline int64_t laplace_idx(int64_t i, int64_t j, int64_t k, int64_t nj, int64_t nk) { return j + nj * ( k + nk * i); }

/*
*  placeholder
* */
void laplace_kernel(int64_t n, const double* x, double* y);

/*
* Builds the algebraic grid node coordinates from the streamwise station
* distribution and wall contours. Cross-stream coordinates use tanh
* clustering from the symmetry plane (s=0) to the wall (s=1).
*
* Flat layout: index = j + nj*(k + nk*i) — j fastest. A caller indexing
* x/y/z directly must declare double buf[ni][nk][nj] to match.
*
* @param ni, nj, nk   Grid dimensions (streamwise, vertical, spanwise)
* @param xs           Streamwise station coordinates, strictly increasing (ni)
* @param hs           Wall half-height at each station (ni)
* @param bs           Wall half-width at each station (ni)
* @param beta_j       Tanh clustering strength, vertical direction (> 0)
* @param beta_k       Tanh clustering strength, spanwise direction (> 0)
* @param x, y, z      Node coordinates, caller-allocated (nj * nk * ni)
* @param ierr         0 = OK; -1 = invalid beta; -2 = non-monotonic xs
*/
void laplace_grid(int64_t ni, int64_t nj, int64_t nk, const double* xs, const double* hs, const double* bs, const double beta_j, const double beta_k, double* x, double* y, double* z, int32_t* ierr);

void laplace_metrics_xi(int64_t ni, int64_t nj, int64_t nk, const double* x, const double* y, const double* z, const double* rho, double* Aii, double* Aij, double* Aik, int32_t* ierr);

void laplace_metrics_eta(int64_t ni, int64_t nj, int64_t nk, const double* x, const double* y, const double* z, const double* rho, double* Aji, double* Ajj, double* Ajk, int32_t* ierr);

void laplace_metrics_zeta(int64_t ni, int64_t nj, int64_t nk, const double* x, const double* y, const double* z, const double* rho, double* Aki, double* Akj, double* Akk, int32_t* ierr);

void laplace_residual(int64_t ni, int64_t nj, int64_t nk, const double* phi, const double* Aii, const double* Aij, const double* Aik, const double* Aji, const double* Ajj, const double* Ajk, const double* Aki, const double* Akj, const double* Akk, double* r, int32_t* ierr);

/*
* Computes the nondimensional density rho/rho_0 at each face (xi, eta, and
* zeta) from the current phi and the grid, via the physical velocity and
* isentropic closure.
*
* @param ni, nj, nk           Grid dimensions (streamwise, vertical, spanwise)
*
* @param x, y, z              Node coordinates (nj, nk, ni)
*
* @param phi                  Potential field (nj, nk, ni)
*
* @param gamma                Ratio of specific heats (> 1)
*
* @param rho_xi               Density at each xi-face (nj, nk, ni-1)
*
* @param rho_eta              Density at each eta-face (nj-1, nk, ni)
*
* @param rho_zeta             Density at each zeta-face (nj, nk-1, ni)
*
*  @param ierr                 0 = OK; -1 = degenerate grid; -2 = sonic limit
*                             reached; -3 = invalid gamma (<= 1)
*/
void laplace_density(int64_t ni, int64_t nj, int64_t nk, const double* x, const double* y, const double* z, const double* phi, const double gamma, double* rho_xi, double* rho_eta, double* rho_zeta, int32_t* ierr);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // LAPLACE_H_

