/* Dependency-free: stdlib + compiler builtins only.
 *
 * Used by laplace_core_solve to publish progress and observe an early-stop
 * request across threads while staying a PURE procedure: a plain
 * intent(inout)/intent(in) Fortran argument gives no guarantee the compiler
 * won't cache a read or reorder a write relative to other memory operations,
 * and Fortran's VOLATILE attribute is explicitly disallowed inside a PURE
 * procedure -- a real atomic release/acquire pair, reached through a
 * pure-asserted bind(C) interface, is the only mechanism that is both
 * standards-legal inside PURE and actually safe across threads. */
#include <stdint.h>

#if defined(__GNUC__) || defined(__clang__)
    #define LAPLACE_ATOMICS_GNU  1
    #define LAPLACE_ATOMICS_MSVC 0
#elif defined(_MSC_VER)
    #define LAPLACE_ATOMICS_GNU  0
    #define LAPLACE_ATOMICS_MSVC 1
#else
    #error "laplace_atomic: unsupported compiler (need GCC/Clang __atomic builtins or MSVC intrinsics)"
#endif

#if LAPLACE_ATOMICS_MSVC
    #include <intrin.h>
    #if defined(_M_X64) || defined(_M_AMD64)
        #define LAPLACE_ATOMICS_BARRIER() _ReadWriteBarrier()
    #elif defined(_M_ARM64)
        #define LAPLACE_ATOMICS_BARRIER() __dmb(0x0B)
    #else
        #error "laplace_atomic: unsupported MSVC architecture (need x64 or ARM64)"
    #endif
#endif

void laplace_atomic_store_rel_i64(int64_t* p, int64_t v)
{
#if LAPLACE_ATOMICS_GNU
    __atomic_store_n(p, v, __ATOMIC_RELEASE);
#else
    LAPLACE_ATOMICS_BARRIER();
    __iso_volatile_store64((volatile __int64*)p, (__int64)v);
#endif
}

/* witness is unused -- its only purpose is to change on every call (the
 * caller passes its loop counter) so a PURE-licensed Fortran compiler can
 * never treat two calls as having "the same arguments" and cache/hoist the
 * load, regardless of what it assumes about p's own constancy. */
int64_t laplace_atomic_load_acq_i64(const int64_t* p, int64_t witness)
{
    (void)witness;
#if LAPLACE_ATOMICS_GNU
    return __atomic_load_n(p, __ATOMIC_ACQUIRE);
#else
    int64_t v = (int64_t)__iso_volatile_load64((const volatile __int64*)p);
    LAPLACE_ATOMICS_BARRIER();
    return v;
#endif
}
