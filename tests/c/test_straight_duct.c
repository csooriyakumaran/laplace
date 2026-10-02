/* Validates README Sec9 #1/#4: for a straight (constant-area) duct, phi = u*x
 * with uniform rho is the EXACT solution -- the discrete residual should be
 * zero to near machine precision, independent of grid stretching, since the
 * metric rule (README Sec4) guarantees exact flux telescoping for uniform flow.
 *
 * Two cases, both exercising the full isentropic closure (not a hardcoded
 * rho=1 shortcut):
 *   - "low-speed"  ~50 m/s at ISA sea-level stagnation conditions (M ~ 0.15)
 *   - "M0.8"       Mach 0.8 directly
 */
#define AETHER_IMPLEMENTATION
#include "aether/aether.h"

#include "laplace/laplace.h"

#include <math.h>
#include <stdio.h>

#define GAMMA 1.4
#define R_AIR 287.05
#define T0_K  288.15

static f64 stagnation_sound_speed(void)
{
    return sqrt(GAMMA * R_AIR * T0_K);
}

/* Isentropic uniform-flow state at Mach M: nondimensional velocity u = phi_x
 * and nondimensional density rho = rho/rho_0 (README Sec1, Sec2). */
static void uniform_flow_state(f64 M, f64* u, f64* rho)
{
    f64 a2 = 1.0 / (1.0 + 0.5 * (GAMMA - 1.0) * M * M);
    *u   = M * sqrt(a2);
    *rho = pow(a2, 1.0 / (GAMMA - 1.0));
}

static int run_case(Arena* arena, const char* label, f64 M,
                     i64 ni, i64 nj, i64 nk,
                     const f64* xs, const f64* hs, const f64* bs,
                     f64 beta_j, f64 beta_k)
{
    i32 ierr = 0;

    f64* x = arena_push_array(arena, f64, nj * nk * ni);
    f64* y = arena_push_array(arena, f64, nj * nk * ni);
    f64* z = arena_push_array(arena, f64, nj * nk * ni);

    laplace_grid(ni, nj, nk, xs, hs, bs, beta_j, beta_k, x, y, z, &ierr);
    if (ierr != 0) { fprintf(stderr, "%s: laplace_grid failed, ierr=%d\n", label, ierr); return 1; }

    f64 u, rho;
    uniform_flow_state(M, &u, &rho);

    f64* phi = arena_push_array(arena, f64, nj * nk * ni);
    for (i64 n = 0; n < nj * nk * ni; ++n) phi[n] = u * x[n];

    f64* rho_xi   = arena_push_array(arena, f64, nj * nk * (ni - 1));
    f64* rho_eta  = arena_push_array(arena, f64, (nj - 1) * nk * ni);
    f64* rho_zeta = arena_push_array(arena, f64, nj * (nk - 1) * ni);
    for (i64 n = 0; n < nj * nk * (ni - 1); ++n) rho_xi[n]   = rho;
    for (i64 n = 0; n < (nj - 1) * nk * ni; ++n) rho_eta[n]  = rho;
    for (i64 n = 0; n < nj * (nk - 1) * ni; ++n) rho_zeta[n] = rho;

    f64* Aii = arena_push_array(arena, f64, nj * nk * (ni - 1));
    f64* Aij = arena_push_array(arena, f64, nj * nk * (ni - 1));
    f64* Aik = arena_push_array(arena, f64, nj * nk * (ni - 1));
    laplace_metrics_xi(ni, nj, nk, x, y, z, rho_xi, Aii, Aij, Aik, &ierr);
    if (ierr != 0) { fprintf(stderr, "%s: laplace_metrics_xi failed, ierr=%d\n", label, ierr); return 1; }

    f64* Aji = arena_push_array(arena, f64, (nj - 1) * nk * ni);
    f64* Ajj = arena_push_array(arena, f64, (nj - 1) * nk * ni);
    f64* Ajk = arena_push_array(arena, f64, (nj - 1) * nk * ni);
    laplace_metrics_eta(ni, nj, nk, x, y, z, rho_eta, Aji, Ajj, Ajk, &ierr);
    if (ierr != 0) { fprintf(stderr, "%s: laplace_metrics_eta failed, ierr=%d\n", label, ierr); return 1; }

    f64* Aki = arena_push_array(arena, f64, nj * (nk - 1) * ni);
    f64* Akj = arena_push_array(arena, f64, nj * (nk - 1) * ni);
    f64* Akk = arena_push_array(arena, f64, nj * (nk - 1) * ni);
    laplace_metrics_zeta(ni, nj, nk, x, y, z, rho_zeta, Aki, Akj, Akk, &ierr);
    if (ierr != 0) { fprintf(stderr, "%s: laplace_metrics_zeta failed, ierr=%d\n", label, ierr); return 1; }

    f64 m_in = rho * u; /* matches README Sec2's m'' = rho*u/(rho_0*a_0); this
                         * is the exact mass flux the uniform solution itself
                         * carries, so the prescribed inlet condition is
                         * consistent with phi = u*x rather than fighting it */

    f64* r = arena_push_array(arena, f64, nj * nk * (ni - 1));
    laplace_residual(ni, nj, nk, phi, y, z, m_in, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, &ierr);
    if (ierr != 0) { fprintf(stderr, "%s: laplace_residual failed, ierr=%d\n", label, ierr); return 1; }

    f64 max_r = 0.0;
    i64 max_idx = 0;
    for (i64 n = 0; n < nj * nk * (ni - 1); ++n) {
        f64 a = fabs(r[n]);
        if (a > max_r) { max_r = a; max_idx = n; }
    }

    const f64 tol = 1e-9;
    b8 pass = max_r < tol;

    printf("%-10s M=%.4f  u=%.6f  rho=%.6f  max|R|=%.3e  (idx %lld)  %s\n",
           label, M, u, rho, max_r, (long long)max_idx, pass ? "PASS" : "FAIL");

    return pass ? 0 : 1;
}

int main(void)
{
    Arena* arena = arena_alloc(MB(16));

    const i64 ni = 11, nj = 9, nk = 7;
    const f64 L = 5.0, H = 1.0, B = 0.6;
    const f64 beta_j = 2.0, beta_k = 1.8;

    f64* xs = arena_push_array(arena, f64, ni);
    f64* hs = arena_push_array(arena, f64, ni);
    f64* bs = arena_push_array(arena, f64, ni);
    for (i64 i = 0; i < ni; ++i) {
        xs[i] = L * (f64)i / (f64)(ni - 1);
        hs[i] = H;
        bs[i] = B;
    }

    f64 M_low = 50.0 / stagnation_sound_speed();

    int status = 0;
    status |= run_case(arena, "low-speed", M_low, ni, nj, nk, xs, hs, bs, beta_j, beta_k);
    status |= run_case(arena, "M0.8",      0.8,   ni, nj, nk, xs, hs, bs, beta_j, beta_k);

    arena_release(arena);

    return status;
}
