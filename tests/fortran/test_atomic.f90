!> Verifies the atomic store/load shim (src/laplace_atomic.c) that
!! laplace_core_solve's progress-reporting/early-stop machinery depends on.
!!
!! Two independent checks:
!!   1. Round-trip correctness: store a value, load it back, confirm it
!!      matches -- catches ABI/linkage/type-size mistakes in the shim itself.
!!   2. Hoisting resistance: a PURE Fortran subroutine calls a PURE-
!!      interfaced external probe (same calling shape as
!!      laplace_atomic_load_acq_i64, including the witness argument) inside
!!      a loop; the probe increments and returns a counter, ignoring its
!!      arguments entirely. If gfortran ever started treating repeated calls
!!      as "same arguments -> same result" and cached/hoisted them, every
!!      iteration would see the same value instead of a strictly increasing
!!      one. This is the actual compiler-behavior guarantee the solve loop's
!!      stop_flag check depends on -- a regression here would otherwise be a
!!      silent, hard-to-diagnose bug, not a compile error.
program test_atomic

    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    use :: laplace_types, only : ik
    use :: laplace_core,  only : laplace_atomic_store_rel_i64, laplace_atomic_load_acq_i64

    implicit none

    interface
        pure function laplace_test_probe_counter(p, witness) result(v) &
                bind(C, name="laplace_test_probe_counter")
            import :: ik
            integer(ik), intent(in)        :: p
            integer(ik), intent(in), value :: witness
            integer(ik)                    :: v
        end function laplace_test_probe_counter
    end interface

    integer(ik), parameter :: n = 200_ik

    integer(ik) :: shared_val
    integer(ik) :: seen(n)
    integer(ik) :: i
    integer     :: status, check_fail

    status = 0

    ! --- Check 1: round-trip correctness -------------------------------------
    check_fail = 0

    shared_val = 0_ik
    call laplace_atomic_store_rel_i64(shared_val, 123456789_ik)
    if (laplace_atomic_load_acq_i64(shared_val, 0_ik) /= 123456789_ik) check_fail = check_fail + 1

    call laplace_atomic_store_rel_i64(shared_val, -42_ik)
    if (laplace_atomic_load_acq_i64(shared_val, 0_ik) /= -42_ik) check_fail = check_fail + 1

    write(output_unit, '(A, A)') 'round-trip       ', merge('PASS', 'FAIL', check_fail == 0)
    status = status + check_fail

    ! --- Check 2: hoisting resistance -----------------------------------------
    check_fail = 0

    call probe_loop(n, 0_ik, seen)

    do i = 1, n
        if (seen(i) /= i) then
            write(error_unit, '(A,I0,A,I0,A,I0)') &
                'hoisting check mismatch at i=', i, '  expected=', i, '  got=', seen(i)
            check_fail = check_fail + 1
        end if
    end do

    write(output_unit, '(A, A)') 'hoisting-resist  ', merge('PASS', 'FAIL', check_fail == 0)
    status = status + check_fail

    if (status /= 0) then
        write(error_unit, '(A)') 'FAILED'
        stop 1
    end if

contains

    pure subroutine probe_loop(n, p, seen)
        integer(ik), intent(in)  :: n
        integer(ik), intent(in)  :: p        ! unchanging per this subroutine's own view
        integer(ik), intent(out) :: seen(n)
        integer(ik) :: i
        do i = 1, n
            seen(i) = laplace_test_probe_counter(p, i)
        end do
    end subroutine probe_loop

end program test_atomic
