!> Validates README Sec9 #1/#4: for a straight (constant-area) duct, phi = u*x
!! with uniform rho is the EXACT solution -- the discrete residual should be
!! zero to near machine precision, independent of grid stretching, since the
!! metric rule (README Sec4) guarantees exact flux telescoping for uniform flow.
!!
!! Two cases, both exercising the full isentropic closure (not a hardcoded
!! rho=1 shortcut):
!!   - "low-speed"  ~50 m/s at ISA sea-level stagnation conditions (M ~ 0.15)
!!   - "M0.8"       Mach 0.8 directly
!!
!! Mirrors tests/c/test_straight_duct.c, but goes through laplace_f_api
!! instead of the C ABI.
program test_straight_duct

    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    use :: laplace_types, only : rk, ik
    use :: laplace_f_api, only : laplace_grid, laplace_metrics_xi, laplace_metrics_eta, &
                                  laplace_metrics_zeta, laplace_residual

    implicit none

    real(rk), parameter :: gamma_ = 1.4_rk
    real(rk), parameter :: r_air  = 287.05_rk
    real(rk), parameter :: t0_k   = 288.15_rk

    integer(ik), parameter :: ni = 11_ik, nj = 9_ik, nk = 7_ik
    real(rk),    parameter :: length = 5.0_rk, half_h = 1.0_rk, half_b = 0.6_rk
    real(rk),    parameter :: beta_j = 2.0_rk, beta_k = 1.8_rk

    real(rk) :: xs(ni), hs(ni), bs(ni)
    real(rk) :: a0, m_low
    integer(ik) :: i
    integer :: status

    do i = 1, ni
        xs(i) = length * real(i - 1, rk) / real(ni - 1, rk)
        hs(i) = half_h
        bs(i) = half_b
    end do

    a0    = sqrt(gamma_ * r_air * t0_k)
    m_low = 50.0_rk / a0

    status = 0
    status = status + run_case('low-speed', m_low,  ni, nj, nk, xs, hs, bs, beta_j, beta_k)
    status = status + run_case('M0.8',      0.8_rk, ni, nj, nk, xs, hs, bs, beta_j, beta_k)

    if (status /= 0) then
        write(error_unit, '(A)') 'FAILED'
        stop 1
    end if

contains

    !> Isentropic uniform-flow state at Mach M: nondimensional velocity
    !! u = phi_xi and density rho = rho/rho_0 (README Sec1, Sec2).
    subroutine uniform_flow_state(m, u, rho)
        real(rk), intent(in)  :: m
        real(rk), intent(out) :: u, rho
        real(rk) :: a2
        a2  = 1.0_rk / (1.0_rk + 0.5_rk * (gamma_ - 1.0_rk) * m * m)
        u   = m * sqrt(a2)
        rho = a2 ** (1.0_rk / (gamma_ - 1.0_rk))
    end subroutine uniform_flow_state

    function run_case(label, m, ni, nj, nk, xs, hs, bs, beta_j, beta_k) result(stat)
        character(*), intent(in) :: label
        real(rk),     intent(in) :: m
        integer(ik),  intent(in) :: ni, nj, nk
        real(rk),     intent(in) :: xs(ni), hs(ni), bs(ni)
        real(rk),     intent(in) :: beta_j, beta_k
        integer :: stat

        real(rk), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:), phi(:,:,:)
        real(rk), allocatable :: rho_xi(:,:,:), rho_eta(:,:,:), rho_zeta(:,:,:)
        real(rk), allocatable :: Aii(:,:,:), Aij(:,:,:), Aik(:,:,:)
        real(rk), allocatable :: Aji(:,:,:), Ajj(:,:,:), Ajk(:,:,:)
        real(rk), allocatable :: Aki(:,:,:), Akj(:,:,:), Akk(:,:,:)
        real(rk), allocatable :: r(:,:,:)

        real(rk) :: u, rho, max_r, tol
        integer  :: ierr
        logical  :: pass

        allocate(x(nj,nk,ni), y(nj,nk,ni), z(nj,nk,ni), phi(nj,nk,ni))

        call laplace_grid(ni, nj, nk, xs, hs, bs, beta_j, beta_k, x, y, z, ierr)
        if (ierr /= 0) then
            write(error_unit, '(A,A,I0)') trim(label), ': laplace_grid failed, ierr=', ierr
            stat = 1
            return
        end if

        call uniform_flow_state(m, u, rho)
        phi = u * x

        allocate(rho_xi(nj, nk, ni-1), rho_eta(nj-1, nk, ni), rho_zeta(nj, nk-1, ni))
        rho_xi   = rho
        rho_eta  = rho
        rho_zeta = rho

        allocate(Aii(nj,nk,ni-1), Aij(nj,nk,ni-1), Aik(nj,nk,ni-1))
        call laplace_metrics_xi(ni, nj, nk, x, y, z, rho_xi, Aii, Aij, Aik, ierr)
        if (ierr /= 0) then
            write(error_unit, '(A,A,I0)') trim(label), ': laplace_metrics_xi failed, ierr=', ierr
            stat = 1
            return
        end if

        allocate(Aji(nj-1,nk,ni), Ajj(nj-1,nk,ni), Ajk(nj-1,nk,ni))
        call laplace_metrics_eta(ni, nj, nk, x, y, z, rho_eta, Aji, Ajj, Ajk, ierr)
        if (ierr /= 0) then
            write(error_unit, '(A,A,I0)') trim(label), ': laplace_metrics_eta failed, ierr=', ierr
            stat = 1
            return
        end if

        allocate(Aki(nj,nk-1,ni), Akj(nj,nk-1,ni), Akk(nj,nk-1,ni))
        call laplace_metrics_zeta(ni, nj, nk, x, y, z, rho_zeta, Aki, Akj, Akk, ierr)
        if (ierr /= 0) then
            write(error_unit, '(A,A,I0)') trim(label), ': laplace_metrics_zeta failed, ierr=', ierr
            stat = 1
            return
        end if

        ! matches README Sec2's m'' = rho*u/(rho_0*a_0); this is the exact mass
        ! flux the uniform solution itself carries, so the prescribed inlet
        ! condition is consistent with phi = u*x rather than fighting it
        allocate(r(nj, nk, ni-1))
        call laplace_residual(ni, nj, nk, phi, y, z, rho * u, &
                               Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, ierr)
        if (ierr /= 0) then
            write(error_unit, '(A,A,I0)') trim(label), ': laplace_residual failed, ierr=', ierr
            stat = 1
            return
        end if

        max_r = maxval(abs(r))
        tol   = 1.0e-9_rk
        pass  = max_r < tol

        write(output_unit, '(A10, A, F8.4, A, F10.6, A, F10.6, A, ES10.3, A, A)') &
            trim(label), '  M=', m, '  u=', u, '  rho=', rho, '  max|R|=', max_r, '  ', &
            merge('PASS', 'FAIL', pass)

        stat = merge(0, 1, pass)

    end function run_case

end program test_straight_duct
