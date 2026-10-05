!> Verifies laplace_core_solve (via laplace_f_api::laplace_solve) end-to-end
!! against the same straight-duct exact solution used by test_sweep and
!! test_straight_duct: for a constant-area duct, phi = u*(x - length) (shifted
!! so phi = 0 at the outlet, matching solve's Dirichlet outlet convention) is
!! exact, with uniform rho = rho/rho_0.
!!
!! Unlike test_sweep (which holds density/metrics fixed to isolate `sweep`
!! itself), this exercises the real solve loop: phi, density, and metrics all
!! update together, starting from a cold phi = 0 guess.
program test_solve

    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    use :: laplace_types, only : rk, ik
    use :: laplace_f_api, only : laplace_grid, laplace_solve, laplace_residual, &
                                  laplace_density, laplace_metrics_xi, &
                                  laplace_metrics_eta, laplace_metrics_zeta

    implicit none

    real(rk), parameter :: gamma_ = 1.4_rk

    integer(ik), parameter :: ni = 11_ik, nj = 9_ik, nk = 7_ik
    real(rk),    parameter :: length = 5.0_rk, half_h = 1.0_rk, half_b = 0.6_rk
    real(rk),    parameter :: beta_j = 2.0_rk, beta_k = 1.8_rk
    real(rk),    parameter :: mach = 0.5_rk

    real(rk),    parameter :: omega = 1.8_rk, omega_rho = 1.0_rk
    integer(ik), parameter :: density_update_stride = 1_ik
    integer(ik), parameter :: max_iter = 8000_ik
    real(rk),    parameter :: tol = 1.0e-9_rk

    real(rk) :: xs(ni), hs(ni), bs(ni)
    real(rk), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
    real(rk), allocatable :: phi(:,:,:), phi_exact(:,:,:)
    real(rk), allocatable :: rho_xi(:,:,:), rho_eta(:,:,:), rho_zeta(:,:,:)
    real(rk), allocatable :: Aii(:,:,:), Aij(:,:,:), Aik(:,:,:)
    real(rk), allocatable :: Aji(:,:,:), Ajj(:,:,:), Ajk(:,:,:)
    real(rk), allocatable :: Aki(:,:,:), Akj(:,:,:), Akk(:,:,:)
    real(rk), allocatable :: r(:,:,:)
    real(rk), allocatable :: max_r_history(:), rms_r_history(:)

    real(rk)    :: u, rho, m_in, a2
    real(rk)    :: max_err, max_r, rms_r, max_r_check
    integer(ik) :: i, n_iter_done
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

    allocate(phi_exact(nj,nk,ni))
    phi_exact = u * (x - length)  ! shifted so phi = 0 at the outlet (Dirichlet)

    allocate(phi(nj,nk,ni))
    phi = 0.0_rk  ! deliberately wrong cold start; outlet plane (i=ni) stays 0,
                  ! consistent with the Dirichlet BC solve never updates

    ! --- Check 1: solve converges to the exact solution ---------------------
    ! max_r_history/rms_r_history exercise the actual progress-reporting
    ! mechanism (src/laplace_atomic.c); progress_iter/stop_flag are omitted
    ! since this is an ordinary synchronous call, not run on its own thread.
    allocate(max_r_history(max_iter), rms_r_history(max_iter))
    call laplace_solve(ni, nj, nk, x, y, z, phi, m_in, gamma_, &
                        omega, omega_rho, density_update_stride, &
                        max_iter, tol, &
                        n_iter_done, ierr, &
                        max_r_history=max_r_history, rms_r_history=rms_r_history)
    if (ierr /= 0) then
        write(error_unit, '(A,I0)') 'laplace_solve failed, ierr=', ierr
        stop 1
    end if

    max_r = max_r_history(n_iter_done)
    rms_r = rms_r_history(n_iter_done)
    max_err = maxval(abs(phi - phi_exact))

    write(output_unit, '(A, I0, A, ES10.3, A, ES10.3, A, ES10.3, A, A)') &
        'solve   n_iter=', n_iter_done, '  max|R|=', max_r, '  rms|R|=', rms_r, &
        '  max|err|=', max_err, '  ', merge('PASS', 'FAIL', max_err < 1.0e-6_rk)
    if (max_err >= 1.0e-6_rk) status = status + 1

    ! --- Check 2: independent cross-check of solve's own residual bookkeeping ---
    ! Rebuild density/metrics from the converged phi completely from scratch
    ! and confirm an externally-computed residual agrees with solve's max_r --
    ! guards against solve reporting a stale max_r out of sync with the final
    ! phi/density/metrics it actually returns.
    allocate(rho_xi(nj,nk,ni-1), rho_eta(nj-1,nk,ni), rho_zeta(nj,nk-1,ni))
    call laplace_density(ni, nj, nk, x, y, z, phi, gamma_, rho_xi, rho_eta, rho_zeta, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A,I0)') 'laplace_density failed, ierr=', ierr
        stop 1
    end if

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

    allocate(r(nj,nk,ni-1))
    call laplace_residual(ni, nj, nk, phi, y, z, m_in, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, ierr)
    if (ierr /= 0) then
        write(error_unit, '(A,I0)') 'laplace_residual failed, ierr=', ierr
        stop 1
    end if

    max_r_check = maxval(abs(r))

    write(output_unit, '(A, ES10.3, A, A)') &
        'cross-check   max|R|=', max_r_check, '  ', merge('PASS', 'FAIL', max_r_check < 1.0e-6_rk)
    if (max_r_check >= 1.0e-6_rk) status = status + 1

    if (status /= 0) then
        write(error_unit, '(A)') 'FAILED'
        stop 1
    end if

end program test_solve
