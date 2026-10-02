#define AETHER_IMPLEMENTATION
#include "aether/aether.h"

#include "laplace/laplace.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv)
{

    if (argc < 3)
    {
        fprintf(stderr, "usage: %s <wall_countour.csv> <grid_out.csv> [nj] [nk] [beta_j] [beta_k]\n", argv[0]);
        return 1;
    }

    i64 nj     = (argc > 3) ? atoll(argv[3]) : 15;
    i64 nk     = (argc > 4) ? atoll(argv[4]) : 15;
    f64 beta_j = (argc > 5) ? atof(argv[5])  : 2.0;
    f64 beta_k = (argc > 6) ? atof(argv[6])  : 2.0;

    Arena* arena = arena_alloc(MB(16));

    bytes file = file_read(arena, argv[1]);
    if (!file.data) { fprintf(stderr, "could not read %s", argv[1]); return 1; }

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

    f64* x = arena_push_array(arena, f64, nj * nk * ni);
    f64* y = arena_push_array(arena, f64, nj * nk * ni);
    f64* z = arena_push_array(arena, f64, nj * nk * ni);

    i32 ierr = 0;

    /* call into fortran kernel code */
    laplace_grid(ni, nj, nk, xs, hs, bs, beta_j, beta_k, x, y, z, &ierr);
    if (ierr != 0) { fprintf(stderr, "laplace_grid faield, ierr=%d\n", ierr); return ierr; }

    /* output grid data */
    FILE* out = fopen(argv[2], "w");
    fprintf(out, "i,j,k,x,y,z\n");
 
    for (i64 i = 0; i < ni; ++i)
    {
        for (i64 k = 0; k < nk; ++k)
        {
            for (i64 j = 0; j < nj; ++j)
            {
                i64 idx = laplace_idx(i, j, k, nj, nk);
                fprintf(out, "%lld,%lld,%lld,%.10f,%.10f,%.10f\n", (long long)i, (long long)j, (long long)k, x[idx], y[idx], z[idx]);
            }
        }
    }

    fclose(out);
    arena_release(arena);

    return 0;
}
