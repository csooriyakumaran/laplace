!> Verifies laplace_core_output (via laplace_f_api::laplace_output) against
!! the straight-duct exact uniform-flow solution used throughout this test
!! suite: for phi = u*x with uniform rho, phi has no j/k dependence at all, so
!! phi_eta = phi_zeta = 0 exactly (not just approximately) and
!! u = phi_xi/x_xi = u_exact exactly, regardless of any truncation in the
!! grid-clustering metric terms themselves -- same cancellation argument as
!! the residual tests' "exact flux telescoping for uniform flow". So
!! velocity/rho/p/T/M should match the closed-form isentropic state to
!! machine precision at every node, boundaries included.
program test_output

    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    use :: laplace_types, only : rk, ik
    use :: laplace_f_api, only : laplace_grid, laplace_output

    implicit none

    real(rk), parameter :: gamma_ = 1.4_rk

    integer(ik), parameter :: ni = 11_ik, nj = 9_ik, nk = 7_ik
    real(rk),    parameter :: length = 5.0_rk, half_h = 1.0_rk, half_b = 0.6_rk
    real(rk),    parameter :: beta_j = 2.0_rk, beta_k = 1.8_rk
    real(rk),    parameter :: mach = 0.5_rk

    real(rk) :: xs(ni), hs(ni), bs(ni)
    real(rk), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:), phi(:,:,:)
    real(rk), allocatable :: u_c(:,:,:), v_c(:,:,:), w_c(:,:,:)
    real(rk), allocatable :: rho_c(:,:,:), p_c(:,:,:), t_c(:,:,:), m_c(:,:,:)

    real(rk)    :: u_x, rho_x, t_x, p_x, a2, max_err
    integer(ik) :: i
    integer     :: ierr, status

    status = 0

    do i = 1, ni
        xs(i) = length * real(i - 1, rk) / real(ni - 1, rk)
        hs(i) = half_h
        bs(i) = half_b
    end do

    allocate(x(nj,nk,ni), y(nj,nk,ni), z(nj,nk,ni))
    call laplace_grid(ni, nj, nk, xs, hs, bs, beta_j, beta_k, x, y, z, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A,I0)') 'laplace_grid failed, ierr=', ierr
        stop 1
    end if

    ! isentropic uniform-flow state at Mach M (same closure as
    ! test_straight_duct.f90/test_solve.f90)
    a2    = 1.0_rk / (1.0_rk + 0.5_rk * (gamma_ - 1.0_rk) * mach * mach)
    u_x   = mach * sqrt(a2)
    rho_x = a2 ** (1.0_rk / (gamma_ - 1.0_rk))
    t_x   = rho_x ** (gamma_ - 1.0_rk)
    p_x   = rho_x ** gamma_

    allocate(phi(nj,nk,ni))
    phi = u_x * x

    allocate(u_c(nj,nk,ni), v_c(nj,nk,ni), w_c(nj,nk,ni))
    allocate(rho_c(nj,nk,ni), p_c(nj,nk,ni), t_c(nj,nk,ni), m_c(nj,nk,ni))

    call laplace_output(ni, nj, nk, x, y, z, phi, gamma_, u_c, v_c, w_c, rho_c, p_c, t_c, m_c, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A,I0)') 'laplace_output failed, ierr=', ierr
        stop 1
    end if

    max_err = max(maxval(abs(u_c - u_x)),   maxval(abs(v_c)),         maxval(abs(w_c)), &
                   maxval(abs(rho_c - rho_x)), maxval(abs(p_c - p_x)), &
                   maxval(abs(t_c - t_x)),     maxval(abs(m_c - mach)))

    write(output_unit, '(A, ES10.3, A, A)') &
        'output   max|err|=', max_err, '  ', merge('PASS', 'FAIL', max_err < 1.0e-9_rk)
    if (max_err >= 1.0e-9_rk) status = status + 1

    if (status /= 0) then
        write(error_unit, '(A)') 'FAILED'
        stop 1
    end if

end program test_output
