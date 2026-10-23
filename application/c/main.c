#define AETHER_IMPLEMENTATION
#include "aether/aether.h"

#include "laplace/laplace.h"

#include "args.h"

#include <math.h>
#include <stdio.h>

#define VERSION_STRING "v0.0.1"

typedef struct Args
{
    f64 M;
    f64 Po;
    f64 To;
    f64 gamma; 
    u32 nj;
    u32 nk;
    f64 beta;
    f64 tolerance;
    b8  help;
    b8  version;
    str8 coords_path;
    str8 output_path;
} Args;

static inline void print_version(void)
{
    printf("Laplace Potential Flow Solver: %s", VERSION_STRING);
}

static inline Args default_args(void)
{
    return (Args){
        .M       = 0.1,
        .To      = 293.15,
        .Po      = 101325.0,
        .gamma   = 1.4,
        .nj      = 15,
        .nk      = 15,
        .beta    = 1.2,
        .tolerance = 1e-6,
        .help    = false,
        .version = false,
        .coords_path = STR("coords.csv"),
        .output_path = STR("out.csv"),
    };
}

static void print_args(Args* args)
{
    fprintf(stdout, "Arguments:\n");
    fprintf(stdout, "------------------------------------\n");
    fprintf(stdout, "  Mach:   %.2f\n", args->M);
    fprintf(stdout, "  Po:     %.2f Pa\n", args->Po);
    fprintf(stdout, "  To:     %.2f K\n", args->To);
    fprintf(stdout, "  gamma:  %.1f\n", args->gamma);
    fprintf(stdout, "  nj:     %lld\n", (long long)args->nj);
    fprintf(stdout, "  nk:     %lld\n", (long long)args->nk);
    fprintf(stdout, "  beta:   %.1f\n", args->beta);
    fprintf(stdout, "  tol:    %.2e\n", args->tolerance);
    fprintf(stdout, "  coords: " STR8_FMT "\n", STR8_ARG(args->coords_path));
    fprintf(stdout, "  output: " STR8_FMT "\n", STR8_ARG(args->output_path));
    fprintf(stdout, "------------------------------------\n");
}

static void print_geom(const f64* xs, const f64* hs, const f64* bs, const i64 ni)
{
    fprintf(stdout, "Contraction Geometry:\n");
    fprintf(stdout, "------------------------------------\n");
    fprintf(stdout, "  ni:     %lld\n",   (long long)ni);
    fprintf(stdout, "  A-in:   %.2f m2\n", hs[0]*bs[0]/(1000*1000));
    fprintf(stdout, "  A-out:  %.2f m2\n", hs[ni-1]*bs[ni-1]/(1000*1000));
    fprintf(stdout, "  CR:     %.2f\n",    hs[0]*bs[0]/hs[ni-1]/bs[ni-1]);
    fprintf(stdout, "  L:      %.2f m\n",  (xs[ni-1] - xs[0])/1000);
    fprintf(stdout, "------------------------------------\n");
}

static b8 load_contour_csv(Arena* scratch, const char* path, f64** out_xs, f64** out_hs, f64** out_bs, i64* out_count, str8* err)
{

    /* read in file data into scratch arena */
    bytes raw = file_read(scratch, path);
    if (!raw.data) { if (err) *err = str8_push_fmt(scratch, "could not read %s", path); return false; }

    /* split file data into lines */
    Str8List lines = str8_split(scratch, (str8){ raw.data, raw.size }, STR("\n"), Str8SplitFlags_SkipEmpty | Str8SplitFlags_Trim);
    i64 ni = (i64)lines.count - 1; /* excluding header */
    if (ni <= 0) { if (err) *err = str8_push_fmt(scratch, "expected a header row plus at least one data row\n"); return false; }

    /* allocate output arrays on scratch arena */
    f64* xs = arena_push_array(scratch, f64, ni);
    f64* hs = arena_push_array(scratch, f64, ni);
    f64* bs = arena_push_array(scratch, f64, ni);

    Str8Node* row = lines.first->next; /* skipping header */

    /* parse lines into output arrays */
    for (i64 i = 0; i < ni; ++i, row = row->next)
    {
        Str8List fields = str8_split(scratch, row->v, STR(","), Str8SplitFlags_Trim);
        if (fields.count != 3) { if (err) *err = str8_push_fmt(scratch, "row %lld: expected 3 columns, got %llu", (long long)i, (unsigned long long)fields.count); return false; }
        Str8Node* f = fields.first;

        if (!str8_to_f64(f->v, &xs[i])) { if (err) *err = str8_push_fmt(scratch, "row %lld: bad x\n", (long long)i); return false; }
        f = f->next;

        if (!str8_to_f64(f->v, &hs[i])) { if (err) *err = str8_push_fmt(scratch, "row %lld: bad y\n", (long long)i); return false; }
        f = f->next;

        if (!str8_to_f64(f->v, &bs[i])) { if (err) *err = str8_push_fmt(scratch, "row %lld: bad x\n", (long long)i); return false; }
    }

    *out_xs = xs; *out_hs = hs; *out_bs = bs; *out_count = ni;
    return true;
}

static f64 generate_scaled_domian(Arena* perm, const f64* xs_c, const f64* hs_c, const f64* bs_c, i64 ni_c, f64** out_xs, f64** out_hs, f64** out_bs, i64* out_ni)
{
    f64 sf     = 1.0 / hs_c[ni_c-1];

    f64 dx_in  = xs_c[1] - xs_c[0];
    f64 dx_out = xs_c[ni_c - 1] - xs_c[ni_c - 2];

    i64 n_up   = (i64)llround(8.0  * hs_c[0]        / dx_in);  /* extend 4 duct width upstream */
    i64 n_down = (i64)llround(16.0 * hs_c[ni_c - 1] / dx_out); /* extend 8 duct width downstream */

    i64 ni     = n_up + ni_c + n_down;

    f64* xs    = arena_push_array(perm, f64, ni);
    f64* hs    = arena_push_array(perm, f64, ni);
    f64* bs    = arena_push_array(perm, f64, ni);

    for (i64 i = 0; i < n_up; ++i)
    {
        xs[i] = sf * xs_c[0] - sf * (f64)(n_up - i) * dx_in ;
        hs[i] = sf * hs_c[0];
        bs[i] = sf * bs_c[0];
    }

    for (i64 i = 0; i < ni_c; ++i)
    {
        xs[n_up + i] = sf * xs_c[i];
        hs[n_up + i] = sf * hs_c[i];
        bs[n_up + i] = sf * bs_c[i];
    }

    for (i64 i = 0; i < n_down; ++i)
    {
        xs[n_up + ni_c + i] = sf * xs_c[ni_c - 1] + sf * (f64)(i + 1) * dx_out;
        hs[n_up + ni_c + i] = sf * hs_c[ni_c - 1];
        bs[n_up + ni_c + i] = sf * bs_c[ni_c - 1];
    }

    *out_xs = xs;
    *out_hs = hs;
    *out_bs = bs;
    *out_ni = ni;
    return sf;
}

/* A/A* = f(M) */
static f64 area_mach_ratio(f64 mach, f64 gamma)
{
    f64 exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0));
    f64 a2       = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;
    return (1.0 / mach) * pow( 2.0 * a2 / (gamma + 1.0), exponent);
}

static f64 subsonic_mach_from_area_ratio(f64 area_ratio, f64 gamma)
{
    f64 lo = 1.0e-5, hi = 1.0 - 1.0e-9;

    for (i32 iter = 0; iter < 60; ++iter)
    {
        f64 mid = 0.5 * (lo + hi);
        if (area_mach_ratio(mid, gamma) > area_ratio) lo = mid;
        else                                     hi = mid;
    }
    return 0.5 * (lo + hi);

}

static inline f64 static_temperature(f64 To, f64 M, f64 k) { return To    /    ( 1.0 + 0.5 * (k - 1.0) * M * M); }
static inline f64 static_pressure(   f64 Po, f64 M, f64 k) { return Po    * pow( 1.0 + 0.5 * (k - 1.0) * M * M, -k   / (k - 1.0) ); }
static inline f64 static_density( f64 rho_o, f64 M, f64 k) { return rho_o * pow( 1.0 + 0.5 * (k - 1.0) * M * M, -1.0 / (k - 1.0) ); }


typedef struct SolveJob
{
    i64 ni, nj, nk;
    const f64 *x, *y, *z;
    f64* phi;
    f64 m_in, gamma, omega, omega_rho, tol;
    i64 density_update_stride, max_iter;
    f64 *max_r_history, *rms_r_history;

    i64 n_iter_done;
    i32 ierr;
    
    u64 progress_iter;
    u64 stop_flag;
    u64 done;

} SolveJob;

static int solve_thread_fn(void* user)
{
    SolveJob* job = (SolveJob*)user;
    laplace_solve(job->ni, job->nj, job->nk, job->x, job->y, job->z, job->phi,
                  job->m_in, job->gamma, job->omega, job->omega_rho,
                  job->density_update_stride, job->max_iter, job->tol,
                  &job->n_iter_done, &job->ierr,
                  (i64*)&job->progress_iter, job->max_r_history, job->rms_r_history,
                  (const i64*)&job->stop_flag);
    atomic_store_rel_u64(&job->done, 1);
    return 0;
}

static void on_interrupt(void* user)
{
    atomic_store_rel_u64(&((SolveJob*)user)->stop_flag, 1);
}

int main(int argc, char** argv)
{
    Args args = default_args();
    {   /* --- Parse Command Line Args --- */
        OptSpec opts[] = {
            { STR("COORDS FILE"), {0},      OptKind_Str,  &args.coords_path, {0},          STR("Input wall contour coordinates CSV      (default: coords.csv)"), true},
            { STR("OUTPUT FILE"), {0},      OptKind_Str,  &args.output_path, {0},          STR("Output CSV                              (default: out.csv)"), true},
            { STR("--help"),    STR("-h"),  OptKind_Flag, &args.help,        {0},          STR("Print help message and exit")},
            { STR("--version"), STR("-v"),  OptKind_Flag, &args.version,     {0},          STR("Print version number and exit")},
            { STR("--Mach"),    STR("-M"),  OptKind_F64,  &args.M,           STR("VALUE"), STR("Test section Mach number                (default: 0.1)")},
            { STR("--Po"),      {0},        OptKind_F64,  &args.Po,          STR("VALUE"), STR("Test section stagnation pressure [Pa]   (default: 101235 Pa)")},
            { STR("--To"),      {0},        OptKind_F64,  &args.To,          STR("VALUE"), STR("Test section stagnation temperature [K] (default: 293.15 K)")},
            { STR("--gamma"),   STR("-k"),  OptKind_F64,  &args.gamma,       STR("VALUE"), STR("Ratio of specific heats                 (default: 1.4)")},
            { STR("--beta"),    STR("-b"),  OptKind_F64,  &args.beta,        STR("VALUE"), STR("Wall-node clustering parameter          (default: 1.2)")},
            { STR("--tol"),     {0},        OptKind_F64,  &args.tolerance,   STR("VALUE"), STR("Convergence tolerance                   (default: 1e-6)")},
            { STR("--nj"),      {0},        OptKind_U32,  &args.nj,          STR("VALUE"), STR("No. of nodes in spanwise direction      (default: 15)")},
            { STR("--nk"),      {0},        OptKind_U32,  &args.nk,          STR("VALUE"), STR("No. of nodes in vertical direction      (default: 15)")},
        };

        if (!parse_args(argc, argv, opts, ARRAY_COUNT(opts))) {                                      return 1; }
        if (args.help)                                        { print_help(opts, ARRAY_COUNT(opts)); return 0; }
        if (args.version)                                     { print_version();                     return 0; }
    }

    print_args(&args);

    /* --- Memory Managemnt --- */
    Arena* memory  = arena_alloc(GB(4));  /* permament memory arena for solution variables */
    Arena* scratch = arena_alloc(GB(4));  /* scratch arena for temporaries -- e.g. when parsing input geometry file, or building the output CSV in memory before a one-shot write */

    /* --- Solver Inputs --- */

    /* flow props */
    const f64 M_ts      = args.M;             /* test section mach number */
    const f64 To        = args.To;            /* stagnation temperature */
    const f64 Po        = args.Po;            /* stagnation pressure */
    const f64 gamma     = args.gamma;         /* ratio of specific heats */
    const f64 R         = 287.05;             /* Gas constant [J/kg-K] */
    /* grid control */
    const i64 nj        = (i64)args.nj;       /* number of nodes in lateral spanwise direction */
    const i64 nk        = (i64)args.nk;       /* number of nodes in vertical spanwise direction */
    const f64 beta      = args.beta;          /* wall-node clustering parameter */
    /* solver parameters */
    const f64 tol       = args.tolerance;
    const f64 omega     = 1.9;                /* SLOR over-relaxation factor */
    const f64 omega_rho = 1.0;                /* Density-update under-relaxation factor */
    const i64 rho_cadence = 10;           /* number of iterations before a density update */
    const i64 max_iter  = 100000;              /* upper bound */

    /* output control */
    const char* coords_filepath = c_str(memory, args.coords_path);
    const char* output_filepath = c_str(memory, args.output_path);
    /* ---------------------------------------------------------------------- */

    f64* xs; /* axial stations */
    f64* hs; /* half-height of input geometry at each axial station */
    f64* bs; /* half-width of input geometry at each axial station */
    i64  ni; /* number of axial stations */
    f64  sf; /* scale factor */

    {   /* --- Domain Generation --- */
        ArenaTemp tmp = arena_begin_temp(scratch);

        f64* xs_c; f64* hs_c; f64* bs_c; i64  ni_c;
        str8 err = {0};

        /* contraction coords loaded from csv: x, h/2, w/2. */
        if (!load_contour_csv(scratch, coords_filepath, &xs_c, &hs_c, &bs_c, &ni_c, &err))
        {
            fprintf(stderr, STR8_FMT, STR8_ARG(err));
            return 1;
        }

        /* report the contraction geometry */
        print_geom(xs_c, hs_c, bs_c, ni_c);

        /* generate the computational domain: extend inlet outlet, normalize on ts width */
        sf = generate_scaled_domian(memory, xs_c, hs_c, bs_c, ni_c, &xs, &hs, &bs, &ni);

        arena_end_temp(tmp);
    }

    /* Grid coordinates: -- memory layout is j + nj * (k + nk * i) */
    f64* x = arena_push_array(memory, f64, nj * nk * ni);
    f64* y = arena_push_array(memory, f64, nj * nk * ni);
    f64* z = arena_push_array(memory, f64, nj * nk * ni);

    { /* --- Computational Grid --- */
        i32 ierr = 0;
        laplace_grid(ni, nj, nk, xs, hs, bs, beta, beta, x, y, z, &ierr);
        if (ierr != 0) { fprintf(stderr, "laplace_grid faield, ierr=%d\n", ierr); return ierr; }
    }

    /* potential feild */
    f64* phi = arena_push_array(memory, f64, nj * nk * ni);

    { /* --- Initial Conditions --- */

        ArenaTemp tmp = arena_begin_temp(scratch);
        f64* u_ao     = arena_push_array(scratch, f64, ni);
        f64* phi_x = arena_push_array(scratch, f64, ni);

        const f64 A_ts     = hs[ni-1] * bs[ni-1]; /* test section area */
        const f64 ratio_ts = area_mach_ratio(M_ts, gamma);

        /* compute local velocity */
        for (i64 i = 0; i < ni; ++i)
        {
            f64 ratio_i = ratio_ts * ( hs[i] * bs[i] ) / A_ts;
            f64 M_i     = subsonic_mach_from_area_ratio(ratio_i, gamma);
            f64 a2_i    = 1.0 / (1.0 + 0.5 * (gamma - 1.0) * M_i * M_i );
            u_ao[i]     = M_i * sqrt(a2_i);
        }

        phi_x[ni-1] = 0.0; /* Dirichlet outlet convention */

        /* u = d/dx(phi) as first approx. */
        for (i64 i = ni - 2; i >= 0; --i)
            phi_x[i] = phi_x[i+1] - 0.5 * (u_ao[i+1] + u_ao[i]) * (xs[i+1]-xs[i]);

        for (i64 i = 0; i < ni; ++i)
            for (i64 n = 0; n < nj * nk; ++n)
                phi[n + nj * nk * i] = phi_x[i];

        /* u and phi_x no longer needed */
        arena_end_temp(tmp);
    }

    i64 n_iter_done = 0;
    f64* max_r_history = arena_push_array(memory, f64, max_iter);
    f64* rms_r_history = arena_push_array(memory, f64, max_iter);

    // todo(chris): move this to separate thread
    {   /* Solver block */

        /* mass flux in/out */
        const f64 A_ts   = hs[ni-1] * bs[ni-1];
        const f64 A_in   = hs[0]    * bs[0];

        const f64 a2_ts  = 1.0 / (1.0 + 0.5 * (gamma - 1.0) * M_ts * M_ts );
        const f64 u_ts   = M_ts * sqrt(a2_ts); /* u/ao */
        const f64 rho_ts = pow( a2_ts, 1.0 / (gamma - 1.0) );
        const f64 m_in   = rho_ts * u_ts * A_ts / A_in;

        SolveJob job = {0};
        job.ni = ni; job.nj = nj; job.nk = nk;
        job.x = x; job.y = y; job.z = z; job.phi = phi;
        job.m_in = m_in; job.gamma = gamma;
        job.omega = omega; job.omega_rho = omega_rho;
        job.density_update_stride = rho_cadence;
        job.max_iter = max_iter; job.tol = tol;
        job.max_r_history = max_r_history; job.rms_r_history = rms_r_history;


        if (!console_signal_install(on_interrupt, &job))
            fprintf(stderr, "warning: could not install Ctrl+C handler\n");

        Thread t = thread_create(solve_thread_fn, &job);
        if (!t.handle) { fprintf(stderr, "failed to start solver thread\n"); return 1; }

        u64 last_reported = 0;
        while (!atomic_load_acq_u64(&job.done))
        {
            thread_sleep_ms(20);
            u64 cur = atomic_load_acq_u64(&job.progress_iter);
            if (cur - last_reported >= 50)
            {
                last_reported = cur - (cur % 50);
                f64 max_r = max_r_history[last_reported - 1];
                f64 rms_r = rms_r_history[last_reported - 1];
                fprintf(stderr, "laplace_solve: iter=%llu  max|R|=%.3e  rms|R|=%.3e\n",
                        (unsigned long long)last_reported, max_r, rms_r);
            }
        }

        int thread_code = 0;
        thread_join(&t, &thread_code);
        console_signal_uninstall();

        n_iter_done  = job.n_iter_done;
        i32 ierr     = job.ierr;

        if (ierr != 0) { fprintf(stderr, "laplace_solve failed, ierr=%d\n", ierr); return ierr; }

        f64 max_r = (n_iter_done > 0) ? max_r_history[n_iter_done - 1] : 0.0;
        f64 rms_r = (n_iter_done > 0) ? rms_r_history[n_iter_done - 1] : 0.0;

        fprintf(stderr, "laplace_solve: n_iter=%lld  max|R|=%.3e  rms|R|=%.3e  ierr=%d\n",
                (long long)n_iter_done, max_r, rms_r, ierr);

    }

    /* node-centered velocity/rho/p/T/M from the converged phi */
    f64* u   = arena_push_array(memory, f64, nj * nk * ni);
    f64* v   = arena_push_array(memory, f64, nj * nk * ni);
    f64* w   = arena_push_array(memory, f64, nj * nk * ni);
    f64* rho = arena_push_array(memory, f64, nj * nk * ni);
    f64* p   = arena_push_array(memory, f64, nj * nk * ni);
    f64* T   = arena_push_array(memory, f64, nj * nk * ni);
    f64* M   = arena_push_array(memory, f64, nj * nk * ni);

    {   /* Report Results */
        i32 ierr = 0;
        laplace_output(ni, nj, nk, x, y, z, phi, gamma, u, v, w, rho, p, T, M, &ierr);
        if (ierr != 0) { fprintf(stderr, "laplace_output failed, ierr=%d\n", ierr); return ierr; }

        f64 a0   = sqrt(gamma * R * To);
        f64 rho0 = Po / R / To;

        for (i64 n = 0; n < nj * nk * ni; ++n)
        {
            x[n]   /= sf; y[n] /= sf; z[n] /= sf;
            u[n]   *= a0; v[n] *= a0; w[n] *= a0;
            phi[n] *= a0 / sf; /* phi ~ a_0 * h_ts (README Sec2) */
            rho[n] *= rho0;
            p[n]   *= Po;
            T[n]   *= To;
        }
    }

    /* output grid + converged potential + physical output -- streamed row by
     * row through a small fixed buffer, instead of building the whole CSV in
     * memory first */
    fprintf(stdout, "Writing data ...");

    u8 stream_buf[64 * 1024];
    FileStream fs = file_stream_open(output_filepath, stream_buf, sizeof(stream_buf));
    if (!file_stream_valid(&fs))
    {
        fprintf(stderr, "failed to open %s\n", output_filepath);
        return 1;
    }

    file_stream_fmt(&fs, 64, "i,j,k,x,y,z,phi,u,v,w,rho,p,t,mach\n");

    for (i64 i = 0; i < ni && file_stream_ok(&fs); ++i)
    {
        for (i64 k = 0; k < nk; ++k)
        {
            for (i64 j = 0; j < nj; ++j)
            {
                i64 idx = laplace_idx(i, j, k, nj, nk);
                file_stream_fmt(&fs, 512,
                    "%lld,%lld,%lld,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
                    (long long)i, (long long)j, (long long)k, x[idx], y[idx], z[idx], phi[idx],
                    u[idx], v[idx], w[idx], rho[idx], p[idx], T[idx], M[idx]);
            }
        }
    }

    if (!file_stream_close(&fs))
    {
        fprintf(stderr, "failed to write %s\n", output_filepath);
        return 1;
    }

    fprintf(stdout, "\n");
    arena_release(memory);
    arena_release(scratch);

    return 0;
}
