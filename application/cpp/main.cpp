#define AETHER_IMPLEMENTATION
#include "aether/aether.h"

#include "laplace/laplace.h"

#include <cstdio>

// Demo driver: proves the C ABI links and runs, with aether owning the
// buffers that cross it. The arena is pushed once, up front; laplace_kernel
// writes into the same buffers in place, so nothing is reallocated per call.
int main(void)
{
    const int64_t n = 5;

    Arena* arena = arena_alloc(MB(1));

    double* x = arena_push_array(arena, double, n);
    double* y = arena_push_array(arena, double, n);

    for (int64_t i = 0; i < n; ++i) {
        x[i] = static_cast<double>(i + 1);
    }

    laplace_kernel(n, x, y);

    std::printf("laplace cpp_app: placeholder kernel smoke test\n");
    for (int64_t i = 0; i < n; ++i) {
        std::printf("%3lld %10.4f %10.4f\n", static_cast<long long>(i), x[i], y[i]);
    }

    arena_release(arena);

    return 0;
}
