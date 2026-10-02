!> Verifies thomas_solve (README Sec7's tridiagonal line-solve primitive) in
!! isolation: build a known solution vector, forward-multiply it through a
!! hand-picked, non-uniform tridiagonal system to get the RHS, then confirm
!! thomas_solve's solve recovers the known vector to near machine precision.
!!
!! Exercises laplace_core directly, not through laplace_f_api -- this is a
!! white-box unit test of a private numerical primitive, not an end-to-end
!! validation of a public feature (c.f. tests/fortran/test_straight_duct.f90).
program test_thomas_solve

    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    use :: laplace_types, only : rk, ik
    use :: laplace_core,  only : thomas_solve

    implicit none

    integer(ik), parameter :: n = 6_ik

    real(rk) :: sub(n), diag(n), sup(n), rhs(n), x_true(n), x(n)
    real(rk) :: max_err, tol
    integer(ik) :: i
    logical :: pass

    ! non-uniform coefficients and a non-trivial solution vector, deliberately
    ! not a simple arithmetic sequence -- avoids a bug that only shows up when
    ! coefficients or the answer happen to be symmetric/uniform
    diag   = [4.0_rk, 5.0_rk, 6.0_rk, 7.0_rk, 5.5_rk, 4.5_rk]
    sub    = [0.0_rk, 1.0_rk, 1.2_rk, 0.9_rk, 1.1_rk, 1.3_rk]
    sup    = [1.5_rk, 1.1_rk, 1.0_rk, 1.3_rk, 0.8_rk, 0.0_rk]
    x_true = [1.0_rk, -2.0_rk, 3.5_rk, -1.5_rk, 2.2_rk, -0.7_rk]

    ! forward-multiply: rhs = A * x_true (A tridiagonal; sub(1), sup(n) unused)
    rhs(1) = diag(1) * x_true(1) + sup(1) * x_true(2)
    do i = 2, n - 1
        rhs(i) = sub(i) * x_true(i - 1) + diag(i) * x_true(i) + sup(i) * x_true(i + 1)
    end do
    rhs(n) = sub(n) * x_true(n - 1) + diag(n) * x_true(n)

    call thomas_solve(n, sub, diag, sup, rhs, x)

    max_err = maxval(abs(x - x_true))
    tol     = 1.0e-10_rk
    pass    = max_err < tol

    write(output_unit, '(A, ES10.3, A, A)') 'thomas_solve  max|error|=', max_err, '  ', &
        merge('PASS', 'FAIL', pass)

    if (.not. pass) then
        write(error_unit, '(A)') 'FAILED'
        stop 1
    end if

end program test_thomas_solve
