#define AETHER_IMPLEMENTATION
#include "aether/aether.h"

#include "laplace/laplace.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define VERSION_STRING "v0.0.1"

typedef enum { OPT_FLAG, OPT_U16, OPT_U32, OPT_F64, OPT_STR } OptKind;

typedef struct OptSpec
{
    str8    name;
    str8    alias;
    OptKind kind;
    void*   out;
    str8    metavar;
    str8    help;
    b8      positional;
} OptSpec;

static inline OptSpec* find_opt(OptSpec* opts, u64 count, str8 key)
{
    if (!opts) return NULL;

    for (u64 i = 0; i < count; ++i)
    {
        if (str8_eq(opts[i].name, key)) return &opts[i];
        if (str8_eq(opts[i].alias, key) && opts[i].alias.size) return &opts[i];
    }

    return NULL;
}

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

static void print_help(OptSpec* opts, u64 count)
{

    printf("\nUSAGE\n");
    printf("  laplace.exe");
    for (u64 i = 0; i < count; ++i)
    {
        const OptSpec* o = &opts[i];
        if (!o->positional) continue;
        printf(" <" STR8_FMT ">", STR8_ARG(o->name));
    }
    for (u64 i = 0; i < count; ++i)
    {
        const OptSpec* o = &opts[i];
         if (o->positional) continue;
        printf(" [");
        printf(STR8_FMT, STR8_ARG(o->name));
        if (o->kind != OPT_FLAG) printf(" <" STR8_FMT ">", STR8_ARG(o->metavar));
        printf("]");
    }
    printf("\n");

    printf("\nARGUMENTS\n");
    for (u64 i = 0; i < count; ++i)
    {
        const OptSpec* o = &opts[i];
        if (!o->positional) continue;
        printf(" " STR8_FMT "\t\t" STR8_FMT " \n", STR8_ARG(o->name), STR8_ARG(o->help));
    }

    printf("\nFLAGS\n");
    for (u64 i = 0; i < count; ++i)
    {
        const OptSpec* o = &opts[i];
        if (o->positional || o->kind != OPT_FLAG) continue;
        if (o->alias.size > 0) printf(" " STR8_FMT ", ", STR8_ARG(o->alias));
        else                   printf("     ");
        printf(STR8_FMT "\t\t" STR8_FMT "\n",
               STR8_ARG(o->name), STR8_ARG(o->help));
    }

    printf("\nOPTIONS\n");
    for (u64 i = 0; i < count; ++i)
    {
        const OptSpec* o = &opts[i];
        if (o->positional || o->kind == OPT_FLAG) continue;
        if (o->alias.size > 0) printf(" " STR8_FMT ", ", STR8_ARG(o->alias));
        else                   printf("     ");
        printf(STR8_FMT " <" STR8_FMT ">\t" STR8_FMT "\n",
               STR8_ARG(o->name), STR8_ARG(o->metavar), STR8_ARG(o->help));
    }
}

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
}

static inline b8 parse_args(int argc, char** argv, Args* a, OptSpec* opts, u64 count)
{

    u64 next_pos = 0;

    for (int i = 1; i < argc; ++i)
    {
        str8 arg = str8_from_c_str(argv[i]);

        /* positional: match the Nth positional-marked OptSpec entry, in table order */
        if (arg.size == 0 || arg.data[0] != '-')
        {
            OptSpec* opt = NULL;
            for (u64 k = 0, seen = 0; k < count; ++k)
            {
                if (!opts[k].positional) continue;
                if (seen == next_pos) { opt = &opts[k]; break; }
                ++seen;
            }
            if (!opt)
            {
                fprintf(stderr, "unexpected positional argument: " STR8_FMT "\n", STR8_ARG(arg));
                return false;
            }
            *(str8*)opt->out = arg;
            ++next_pos;
            continue;
        }

        /* keyword: --name / --name=vale / --name vale, matched by name or alias */
        str8 key, inline_val;
        b8 has_inline = str8_cut(arg, STR("="), &key, &inline_val);

        OptSpec* opt = find_opt(opts, count, key);
        if (!opt) { fprintf(stderr, "Unrecognized option: " STR8_FMT "\n", STR8_ARG(key)); return false; }

        str8 val = {0};

        if (opt->kind != OPT_FLAG)
        {
            if (has_inline) val = inline_val;
            else if (i + 1 < argc) val = str8_from_c_str(argv[++i]);
            else { fprintf(stderr, STR8_FMT " requires a value\n", STR8_ARG(key)); return false; }
        }

        switch (opt->kind)
        {
            case OPT_FLAG: *(b8*)opt->out = true; break;
            case OPT_U16: {
                u64 v;
                if (!str8_to_u64(val, &v) || v == 0 || v > U16_MAX) {
                    fprintf(stderr, STR8_FMT " value must fit in a u16, passed \n" STR8_FMT, STR8_ARG(key), STR8_ARG(val));
                    return false;
                }
                *(u16*)opt->out = (u16)v;
            } break;
            case OPT_U32:{
                u64 v;
                if (!str8_to_u64(val, &v) || v == 0 || v > U32_MAX) {
                    fprintf(stderr, STR8_FMT " value must fit in a u32, passed \n" STR8_FMT, STR8_ARG(key), STR8_ARG(val));
                    return false;
                }
                *(u32*)opt->out = (u32)v;
            } break;
            case OPT_F64: {
                f64 v; 
                if (!str8_to_f64(val, &v))
                {
                    fprintf(stderr, STR8_FMT " value must be a valid f64, passed \n" STR8_FMT, STR8_ARG(key), STR8_ARG(val));
                    return false;
                }
                *(f64*)opt->out = v;
            } break;
            case OPT_STR: {
                *(str8*)opt->out = val;
            } break;

            default: {
                fprintf(stderr, "Unknown option kind\n");
                return false;
            }
        }
    }

    return true;
}

int main(int argc, char** argv)
{
    Args args = default_args();
    OptSpec opts[] = {
        { STR("COORDS FILE"), {0},      OPT_STR,  &args.coords_path, {0},        STR("Input wall contour coordinates CSV      (default: coords.csv)"), true},
        { STR("OUTPUT FILE"), {0},      OPT_STR,  &args.output_path, {0},        STR("Output CSV                              (default: out.csv)"), true},
        { STR("--help"),    STR("-h"),  OPT_FLAG, &args.help,      {0},          STR("Print help message and exit")},
        { STR("--version"), STR("-v"),  OPT_FLAG, &args.version,   {0},          STR("Print version number and exit")},
        { STR("--Mach"),    STR("-M"),  OPT_F64,  &args.M,         STR("VALUE"), STR("Test section Mach number                (default: 0.1)")},
        { STR("--Po"),      STR("-Po"), OPT_F64,  &args.Po,        STR("VALUE"), STR("Test section stagnation pressure [Pa]   (default: 101235 Pa)")},
        { STR("--To"),      STR("-To"), OPT_F64,  &args.To,        STR("VALUE"), STR("Test section stagnation temperature [K] (default: 293.15 K)")},
        { STR("--gamma"),   STR("-k"),  OPT_F64,  &args.gamma,     STR("VALUE"), STR("Ratio of specific heats                 (default: 1.4)")},
        { STR("--beta"),    STR("-b"),  OPT_F64,  &args.beta,      STR("VALUE"), STR("Wall-node clustering parameter          (default: 1.2)")},
        { STR("--tol"),     {0},        OPT_F64,  &args.tolerance, STR("VALUE"), STR("Convergence tolerance                   (default: 1e-6)")},
        { STR("--nj"),      STR("-nj"), OPT_U32,  &args.nj,        STR("VALUE"), STR("No. of nodes in spanwise direction      (default: 15)")},
        { STR("--nk"),      STR("-nk"), OPT_U32,  &args.nk,        STR("VALUE"), STR("No. of nodes in vertical direction      (default: 15)")},
    };

    if (!parse_args(argc, argv, &args, opts, ARRAY_COUNT(opts))) return 1;
    if (args.help)    { print_help(opts, ARRAY_COUNT(opts));    return 0; }
    if (args.version) { print_version(); return 0; }

    print_args(&args);


    f64 mach   = args.M;
    i64 nj     = (i64)args.nj;
    i64 nk     = (i64)args.nk;
    f64 beta_j = args.beta;
    f64 beta_k = args.beta;
    f64 gamma  = args.gamma;

    Arena* arena = arena_alloc(MB(16));
    const char* coords_filepath = c_str(arena, args.coords_path);
    const char* output_filepath = c_str(arena, args.output_path);

    bytes file = file_read(arena, coords_filepath);
    if (!file.data) { fprintf(stderr, "could not read %s", coords_filepath); return 1; }

    str8 content = {file.data, file.size};
    Str8List lines = str8_split(arena, content, STR("\n"), Str8SplitFlags_SkipEmpty | Str8SplitFlags_Trim);

    if (lines.count < 2) { fprintf(stderr, "expected a header row plus at least one data row\n"); return 1; }

    i64 ni = (i64)lines.count - 1; /* excluding header */

    f64* xs = arena_push_array(arena, f64, ni);
    f64* hs = arena_push_array(arena, f64, ni);
    f64* bs = arena_push_array(arena, f64, ni);

    /* parse the input data */
    Str8Node* row = lines.first->next; /* skip header */
    for (i64 i = 0; i < ni; ++i, row = row->next)
    {
        Str8List fields = str8_split(arena, row->v, STR(","), Str8SplitFlags_Trim);
        if (fields.count != 3)
        {
            fprintf(stderr, "row %lld: expected 3 columns, got %llu", (long long)i, (unsigned long long)fields.count);
            return 1;
        }
        Str8Node* f = fields.first;
        if (!str8_to_f64(f->v, &xs[i])) { fprintf(stderr, "row %lld: bad xs\n", (long long)i); return 1; }
        f = f->next;
    if (!str8_to_f64(f->v, &hs[i])) { fprintf(stderr, "row %lld: bad hs\n", (long long)i); return 1; }
        f = f->next;
        if (!str8_to_f64(f->v, &bs[i])) { fprintf(stderr, "row %lld: bad bs\n", (long long)i); return 1; }
    }

    /* nondimensionalize lengths by the exit (test-section) half-height, per
     * README Sec2 -- the solver's tolerances/omega were tuned assuming
     * O(1)-scale lengths, not raw CSV units (e.g. mm) */
    f64 l_ref = hs[ni - 1];
    for (i64 i = 0; i < ni; ++i) { xs[i] /= l_ref; hs[i] /= l_ref; bs[i] /= l_ref; }
    fprintf(stderr, "solving nondimensionally (l_ref=%.6f, exit half-height); output is rescaled back to physical units\n", l_ref);

    f64* x = arena_push_array(arena, f64, nj * nk * ni);
    f64* y = arena_push_array(arena, f64, nj * nk * ni);
    f64* z = arena_push_array(arena, f64, nj * nk * ni);

    i32 ierr = 0;

    /* call into fortran kernel code */
    laplace_grid(ni, nj, nk, xs, hs, bs, beta_j, beta_k, x, y, z, &ierr);
    if (ierr != 0) { fprintf(stderr, "laplace_grid faield, ierr=%d\n", ierr); return ierr; }

    /* README Sec2: the caller-supplied Mach number is the TARGET TEST-SECTION
     * Mach M_ts (the constant-area end -- index ni-1 here), not the inlet's.
     * m'' at the test section scales to the inlet by the area ratio:
     * m''_in = m'' * A_ts/A_in. Same isentropic closure as
     * tests/fortran/test_solve.f90, evaluated at the test section. */
    f64 a2_ts  = 1.0 / (1.0 + 0.5 * (gamma - 1.0) * mach * mach);
    f64 u_ts   = mach * sqrt(a2_ts);
    f64 rho_ts = pow(a2_ts, 1.0 / (gamma - 1.0));
    f64 mpp_ts = rho_ts * u_ts;

    f64 a_ts = hs[ni - 1] * bs[ni - 1];
    f64 a_in = hs[0]      * bs[0];
    f64 m_in = mpp_ts * a_ts / a_in;

    f64* phi = arena_push_array(arena, f64, nj * nk * ni);
    for (i64 i = 0; i < ni; ++i)
    {
        /* u_ts as a crude constant slope is just a starting guess -- SLOR
         * doesn't need an accurate initial field, only a smooth one */
        f64 phi_i = u_ts * (xs[i] - xs[ni - 1]); /* 0 at the outlet, matching solve's Dirichlet BC */
        for (i64 n = 0; n < nj * nk; ++n) phi[n + nj * nk * i] = phi_i;
    }

    /* SLOR/density-update knobs validated against the straight-duct exact
     * solution in tests/fortran/test_solve.f90; tune max_iter/tol if this
     * contraction needs more (or converges well before) that budget */
    const f64 omega     = 1.5;
    const f64 omega_rho = 1.0;
    const i64 density_update_stride = 10;
    const i64 max_iter  = 8000;
    const f64 tol       = 1.0e-9;

    /* laplace_solve is a pure routine now -- no in-loop progress printing.
     * max_r_history/rms_r_history (1 slot per sweep, caller-allocated) carry
     * what print_every used to; progress_iter/stop_flag (both NULL here)
     * are for a caller running this on its own thread, see laplace.h. */
    i64 n_iter_done = 0;
    f64* max_r_history = arena_push_array(arena, f64, max_iter);
    f64* rms_r_history = arena_push_array(arena, f64, max_iter);

    laplace_solve(ni, nj, nk, x, y, z, phi, m_in, gamma,
                   omega, omega_rho, density_update_stride, max_iter, tol,
                   &n_iter_done, &ierr,
                   NULL, max_r_history, rms_r_history, NULL);

    /* history arrays are 0-indexed in C; sweep n_iter_done's result lives at
     * n_iter_done - 1 */
    f64 max_r = (n_iter_done > 0) ? max_r_history[n_iter_done - 1] : 0.0;
    f64 rms_r = (n_iter_done > 0) ? rms_r_history[n_iter_done - 1] : 0.0;

    fprintf(stderr, "laplace_solve: n_iter=%lld  max|R|=%.3e  rms|R|=%.3e  ierr=%d\n",
            (long long)n_iter_done, max_r, rms_r, ierr);
    if (ierr != 0) { fprintf(stderr, "laplace_solve failed, ierr=%d\n", ierr); return ierr; }

    /* node-centered velocity/rho/p/T/M (README Sec8) from the converged phi */
    f64* u_out = arena_push_array(arena, f64, nj * nk * ni);
    f64* v_out = arena_push_array(arena, f64, nj * nk * ni);
    f64* w_out = arena_push_array(arena, f64, nj * nk * ni);
    f64* rho_out = arena_push_array(arena, f64, nj * nk * ni);
    f64* p_out = arena_push_array(arena, f64, nj * nk * ni);
    f64* t_out = arena_push_array(arena, f64, nj * nk * ni);
    f64* mach_out = arena_push_array(arena, f64, nj * nk * ni);

    laplace_output(ni, nj, nk, x, y, z, phi, gamma, u_out, v_out, w_out, rho_out, p_out, t_out, mach_out, &ierr);
    if (ierr != 0) { fprintf(stderr, "laplace_output failed, ierr=%d\n", ierr); return ierr; }

    /* scale back to physical dimensions (README Sec2) for everything written
     * out -- the solve itself stays nondimensional throughout (that's what
     * tol/omega above are tuned against); only the CSV boundary is physical.
     * Stagnation conditions aren't a CLI input yet, so ISA sea level is
     * assumed here, same constants as tests/c/test_straight_duct.c. */
    const f64 r_air = 287.05;
    const f64 t0_k  = 288.15;
    const f64 p0_pa = 101325.0;
    f64 a0   = sqrt(gamma * r_air * t0_k);
    f64 rho0 = p0_pa / (r_air * t0_k);

    for (i64 n = 0; n < nj * nk * ni; ++n)
    {
        x[n] *= l_ref; y[n] *= l_ref; z[n] *= l_ref;
        phi[n] *= a0 * l_ref; /* phi ~ a_0 * h_ts (README Sec2) */
        u_out[n] *= a0; v_out[n] *= a0; w_out[n] *= a0;
        rho_out[n] *= rho0;
        p_out[n]   *= p0_pa;
        t_out[n]   *= t0_k;
        /* mach_out is already dimensionless -- left as-is */
    }

    /* output grid + converged potential + physical output */
    FILE* out = fopen(output_filepath, "w");
    fprintf(out, "i,j,k,x,y,z,phi,u,v,w,rho,p,t,mach\n");

    for (i64 i = 0; i < ni; ++i)
    {
        for (i64 k = 0; k < nk; ++k)
        {
            for (i64 j = 0; j < nj; ++j)
            {
                i64 idx = laplace_idx(i, j, k, nj, nk);
                fprintf(out, "%lld,%lld,%lld,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
                        (long long)i, (long long)j, (long long)k, x[idx], y[idx], z[idx], phi[idx],
                        u_out[idx], v_out[idx], w_out[idx], rho_out[idx], p_out[idx], t_out[idx], mach_out[idx]);
            }
        }
    }

    fclose(out);
    arena_release(arena);

    return 0;
}
