!> Verifies `sweep` (one SLOR line-relaxation pass, README Sec7) in isolation,
!! before it's wired into the full laplace_core_solve convergence loop.
!!
!! Two checks, both on the straight duct (phi = u*(x - L) is exact, shifted so
!! phi = 0 at the outlet x = L, matching the Dirichlet outlet convention sweep
!! assumes -- phi's own derivatives, and hence rho/metrics, are unaffected by
!! this additive shift):
!!
!!   1. Fixed point: starting from the exact solution, one sweep should leave
!!      phi essentially unchanged (residual's already ~0, so the correction
!!      solved for should be too).
!!   2. Convergence: starting from phi = 0 (a deliberately wrong guess), with
!!      density/metrics held fixed at their correct values, repeated sweeps
!!      should drive phi toward the exact solution.
!!
!! Density/metrics are NOT updated between sweeps here -- that cadence belongs
!! to the outer solve loop, not yet built. This isolates `sweep` itself.
program test_sweep

    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    use :: laplace_types, only : rk, ik
    use :: laplace_core,  only : sweep
    use :: laplace_f_api, only : laplace_grid, laplace_metrics_xi, laplace_metrics_eta, &
                                  laplace_metrics_zeta, laplace_residual

    implicit none

    real(rk), parameter :: gamma_ = 1.4_rk

    integer(ik), parameter :: ni = 11_ik, nj = 9_ik, nk = 7_ik
    real(rk),    parameter :: length = 5.0_rk, half_h = 1.0_rk, half_b = 0.6_rk
    real(rk),    parameter :: beta_j = 2.0_rk, beta_k = 1.8_rk
    real(rk),    parameter :: mach = 0.5_rk

    integer(ik), parameter :: n_iter = 3000_ik ! 500 wasn't enough sweeps to clear the 1e-6 threshold at this omega
    real(rk),    parameter :: omega  = 1.8_rk

    real(rk) :: xs(ni), hs(ni), bs(ni)
    real(rk), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
    real(rk), allocatable :: phi_exact(:,:,:), phi(:,:,:)
    real(rk), allocatable :: rho_xi(:,:,:), rho_eta(:,:,:), rho_zeta(:,:,:)
    real(rk), allocatable :: Aii(:,:,:), Aij(:,:,:), Aik(:,:,:)
    real(rk), allocatable :: Aji(:,:,:), Ajj(:,:,:), Ajk(:,:,:)
    real(rk), allocatable :: Aki(:,:,:), Akj(:,:,:), Akk(:,:,:)
    real(rk), allocatable :: r(:,:,:)

    real(rk)    :: u, rho, m_in, a2
    real(rk)    :: max_move, max_err, max_r
    real(rk)    :: sweep_max_r, sweep_rms_r
    integer(ik) :: i, iter
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

    a2   = 1.0_rk / (1.0_rk + 0.5_rk * (gamma_ - 1.0_rk) * mach * mach)
    u    = mach * sqrt(a2)
    rho  = a2 ** (1.0_rk / (gamma_ - 1.0_rk))
    m_in = rho * u

    allocate(rho_xi(nj,nk,ni-1), rho_eta(nj-1,nk,ni), rho_zeta(nj,nk-1,ni))
    rho_xi = rho; rho_eta = rho; rho_zeta = rho

    allocate(Aii(nj,nk,ni-1), Aij(nj,nk,ni-1), Aik(nj,nk,ni-1))
    call laplace_metrics_xi(ni, nj, nk, x, y, z, rho_xi, Aii, Aij, Aik, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A,I0)') 'laplace_metrics_xi failed, ierr=', ierr
        stop 1
    end if

    allocate(Aji(nj-1,nk,ni), Ajj(nj-1,nk,ni), Ajk(nj-1,nk,ni))
    call laplace_metrics_eta(ni, nj, nk, x, y, z, rho_eta, Aji, Ajj, Ajk, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A,I0)') 'laplace_metrics_eta failed, ierr=', ierr
        stop 1
    end if

    allocate(Aki(nj,nk-1,ni), Akj(nj,nk-1,ni), Akk(nj,nk-1,ni))
    call laplace_metrics_zeta(ni, nj, nk, x, y, z, rho_zeta, Aki, Akj, Akk, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A,I0)') 'laplace_metrics_zeta failed, ierr=', ierr
        stop 1
    end if

    allocate(phi_exact(nj,nk,ni))
    phi_exact = u * (x - length)  ! shifted so phi = 0 at the outlet (Dirichlet)

    ! --- Check 1: fixed point ------------------------------------------------
    allocate(phi(nj,nk,ni))
    phi = phi_exact

    call sweep(ni, nj, nk, phi, y, z, m_in, omega, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, sweep_max_r, sweep_rms_r)

    max_move = maxval(abs(phi - phi_exact))
    write(output_unit, '(A, ES10.3, A, A)') 'fixed-point   max|move|=', max_move, '  ', &
        merge('PASS', 'FAIL', max_move < 1.0e-9_rk)
    if (max_move >= 1.0e-9_rk) status = status + 1

    ! --- Check 2: convergence from a wrong guess -----------------------------
    phi = 0.0_rk

    do iter = 1, n_iter
        call sweep(ni, nj, nk, phi, y, z, m_in, omega, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, sweep_max_r, sweep_rms_r)
    end do

    max_err = maxval(abs(phi - phi_exact))

    allocate(r(nj, nk, ni-1))
    call laplace_residual(ni, nj, nk, phi, y, z, m_in, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A,I0)') 'laplace_residual failed, ierr=', ierr
        stop 1
    end if
    max_r = maxval(abs(r))

    write(output_unit, '(A, ES10.3, A, ES10.3, A, I0, A, A)') &
        'convergence   max|err|=', max_err, '  max|R|=', max_r, '  after ', n_iter, ' sweeps  ', &
        merge('PASS', 'FAIL', max_err < 1.0e-6_rk)
    if (max_err >= 1.0e-6_rk) status = status + 1

    if (status /= 0) then
        write(error_unit, '(A)') 'FAILED'
        stop 1
    end if

end program test_sweep
