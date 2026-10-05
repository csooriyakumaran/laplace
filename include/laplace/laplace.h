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

void laplace_residual(int64_t ni, int64_t nj, int64_t nk, const double* phi, const double* y, const double* z, const double m_in, const double* Aii, const double* Aij, const double* Aik, const double* Aji, const double* Ajj, const double* Ajk, const double* Aki, const double* Akj, const double* Akk, double* r, int32_t* ierr);

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

/*
* Iteratively solves the full potential equation via SLOR line relaxation
* (one `sweep` per iteration) with lagged, under-relaxed density updates,
* starting from the caller's initial guess for phi and converging until
* max|R| < tol or max_iter sweeps have run. All density/metric workspace is
* allocated internally -- the caller only owns the grid and phi.
*
* The last four parameters are all optional -- pass NULL for any/all of them
* for an ordinary synchronous solve. They exist for a caller running this on
* its own thread:
*
*   - progress_iter, if non-NULL, is updated via an atomic release-store
*     every iteration. Poll it from another thread via an acquire-load
*     (e.g. aether's atomic_load_acq_u64); observing N there guarantees
*     max_r_history[0..N-1]/rms_r_history[0..N-1] are safe to read, even if
*     the solve has since raced ahead further.
*   - max_r_history/rms_r_history, if non-NULL, must be caller-allocated
*     arrays of length max_iter; slot i holds the residual norms after
*     iteration i+1 (one iteration behind the just-applied correction,
*     except the final slot, which is corrected to be exact).
*   - stop_flag, if non-NULL, is read every iteration via an atomic
*     acquire-load. Set it nonzero via a release-store to request early
*     exit; phi is left at whatever the last completed sweep produced, same
*     as any other exit.
*
* @param ni, nj, nk              Grid dimensions (streamwise, vertical, spanwise)
* @param x, y, z                 Node coordinates (nj, nk, ni)
* @param phi                     Potential field (nj, nk, ni); initial guess
*                                 in, converged field out
* @param m_in                    Prescribed nondimensional inlet mass flux
* @param gamma                   Ratio of specific heats (> 1)
* @param omega                   SLOR over-relaxation factor
* @param omega_rho               Density-update under-relaxation factor
* @param density_update_stride   Recompute density/metrics every N sweeps
* @param max_iter                Maximum number of sweeps
* @param tol                     Convergence tolerance on max|R|
* @param n_iter_done             Number of sweeps actually performed
* @param ierr                    0 = OK; -1 = degenerate grid; -2 = sonic
*                                 limit reached; -3 = invalid gamma (<= 1)
* @param progress_iter           Optional (NULL-able); see above
* @param max_r_history           Optional (NULL-able), length max_iter; see above
* @param rms_r_history           Optional (NULL-able), length max_iter; see above
* @param stop_flag               Optional (NULL-able); see above
*/
void laplace_solve(int64_t ni, int64_t nj, int64_t nk, const double* x, const double* y, const double* z, double* phi, const double m_in, const double gamma, const double omega, const double omega_rho, const int64_t density_update_stride, const int64_t max_iter, const double tol, int64_t* n_iter_done, int32_t* ierr, int64_t* progress_iter, double* max_r_history, double* rms_r_history, const int64_t* stop_flag);

/*
* Node-centered physical output (README Sec8) at every node of the full
* domain: velocity from node-centered (one-sided at boundaries) derivatives
* of phi, then rho/p/T/M from the same isentropic closure as laplace_density.
* No size reduction -- unlike residual/metrics, these are meaningful at every
* node including boundaries. All nondimensional (README Sec1/Sec2);
* dimensionalizing by rho_0/p_0/T_0/a_0 is the caller's job.
*
* @param ni, nj, nk      Grid dimensions (streamwise, vertical, spanwise)
* @param x, y, z         Node coordinates (nj, nk, ni)
* @param phi             Converged potential field (nj, nk, ni)
* @param gamma           Ratio of specific heats (> 1)
* @param u, v, w         Nondimensional physical velocity (nj, nk, ni)
* @param rho, p, t       Nondimensional density, pressure, temperature (nj, nk, ni)
* @param mach            Mach number (nj, nk, ni)
* @param ierr            0 = OK; -1 = degenerate grid; -2 = sonic limit
*                         reached; -3 = invalid gamma (<= 1)
*/
void laplace_output(int64_t ni, int64_t nj, int64_t nk, const double* x, const double* y, const double* z, const double* phi, const double gamma, double* u, double* v, double* w, double* rho, double* p, double* t, double* mach, int32_t* ierr);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // LAPLACE_H_

