/* Test-only fixture for tests/fortran/test_atomic.f90 -- not part of
 * laplacelib. Increments a counter on every call, ignoring its arguments, to
 * confirm a PURE-interfaced external call made inside a Fortran loop (the
 * same pattern src/laplace_atomic.c's load function uses) is actually
 * invoked fresh every time, not hoisted/cached by the compiler. */
#include <stdint.h>

static int64_t call_count = 0;

int64_t laplace_test_probe_counter(const int64_t* p, int64_t witness)
{
    (void)p;
    (void)witness;
    call_count += 1;
    return call_count;
}
