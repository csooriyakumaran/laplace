module laplace_f_api

    use :: laplace_types, only : rk, ik
    use :: laplace_core,  only : laplace_core_kernel
    use :: laplace_core,  only : laplace_core_grid
    use :: laplace_core,  only : laplace_core_metrics_xi
    use :: laplace_core,  only : laplace_core_metrics_eta
    use :: laplace_core,  only : laplace_core_metrics_zeta
    use :: laplace_core,  only : laplace_core_residual
    use :: laplace_core,  only : laplace_core_density
    use :: laplace_core,  only : laplace_core_solve
    use :: laplace_core,  only : laplace_core_output

    implicit none

    private

    public :: laplace_kernel
    public :: laplace_grid
    public :: laplace_metrics_xi
    public :: laplace_metrics_eta
    public :: laplace_metrics_zeta
    public :: laplace_residual
    public :: laplace_density
    public :: laplace_solve
    public :: laplace_output

contains

    subroutine laplace_kernel(n, x, y)
        integer(ik), intent(in)                :: n
        real(rk),    intent(in),  dimension(n) :: x
        real(rk),    intent(out), dimension(n) :: y
        call laplace_core_kernel(n, x, y)
    end subroutine laplace_kernel

    pure subroutine laplace_grid(ni, nj, nk, xs, hs, bs, beta_j, beta_k, x, y, z, ierr)
        integer(ik), intent(in)                         :: ni, nj, nk
        real(rk),    intent(in),  dimension(ni)         :: xs, hs, bs
        real(rk),    intent(in)                         :: beta_j, beta_k
        real(rk),    intent(out), dimension(nj, nk, ni) :: x, y, z
        integer,     intent(out)                        :: ierr
        call laplace_core_grid(ni, nj, nk, xs, hs, bs, beta_j, beta_k, x, y, z, ierr)
    end subroutine laplace_grid

    pure subroutine laplace_metrics_xi(ni, nj, nk, x, y, z, rho_xi, Aii, Aij, Aik, ierr)
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk, ni)   :: x, y, z
        real(rk),    intent(in),  dimension(nj, nk, ni-1) :: rho_xi
        real(rk),    intent(out), dimension(nj, nk, ni-1) :: Aii, Aij, Aik
        integer,     intent(out)                          :: ierr
        call laplace_core_metrics_xi(ni, nj, nk, x, y, z, rho_xi, Aii, Aij, Aik, ierr)
    end subroutine laplace_metrics_xi

    pure subroutine laplace_metrics_eta(ni, nj, nk, x, y, z, rho_eta, Aji, Ajj, Ajk, ierr)
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj,   nk, ni) :: x, y, z
        real(rk),    intent(in),  dimension(nj-1, nk, ni) :: rho_eta
        real(rk),    intent(out), dimension(nj-1, nk, ni) :: Aji, Ajj, Ajk
        integer,     intent(out)                          :: ierr
        call laplace_core_metrics_eta(ni, nj, nk, x, y, z, rho_eta, Aji, Ajj, Ajk, ierr)
    end subroutine laplace_metrics_eta

    pure subroutine laplace_metrics_zeta(ni, nj, nk, x, y, z, rho_zeta, Aki, Akj, Akk, ierr)
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk,   ni) :: x, y, z
        real(rk),    intent(in),  dimension(nj, nk-1, ni) :: rho_zeta
        real(rk),    intent(out), dimension(nj, nk-1, ni) :: Aki, Akj, Akk
        integer,     intent(out)                          :: ierr
        call laplace_core_metrics_zeta(ni, nj, nk, x, y, z, rho_zeta, Aki, Akj, Akk, ierr)
    end subroutine laplace_metrics_zeta

    pure subroutine laplace_residual(ni, nj, nk, phi, y, z, m_in, &
                                      Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, ierr)
        integer(ik), intent(in)                               :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj,   nk,   ni)   :: phi, y, z
        real(rk),    intent(in)                                :: m_in
        real(rk),    intent(in),  dimension(nj,   nk,   ni-1) :: Aii, Aij, Aik
        real(rk),    intent(in),  dimension(nj-1, nk,   ni)   :: Aji, Ajj, Ajk
        real(rk),    intent(in),  dimension(nj,   nk-1, ni)   :: Aki, Akj, Akk
        real(rk),    intent(out), dimension(nj,   nk,   ni-1) :: r
        integer,     intent(out)                               :: ierr
        call laplace_core_residual(ni, nj, nk, phi, y, z, m_in, &
                                    Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, ierr)
    end subroutine laplace_residual

    pure subroutine laplace_density(ni, nj, nk, x, y, z, phi, gamma, rho_xi, rho_eta, rho_zeta, ierr)
        integer(ik), intent(in)                             :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk, ni)     :: x, y, z, phi
        real(rk),    intent(in)                             :: gamma
        real(rk),    intent(out), dimension(nj,   nk, ni-1) :: rho_xi
        real(rk),    intent(out), dimension(nj-1, nk,   ni) :: rho_eta
        real(rk),    intent(out), dimension(nj,   nk-1, ni) :: rho_zeta
        integer,     intent(out)                            :: ierr
        call laplace_core_density(ni, nj, nk, x, y, z, phi, gamma, rho_xi, rho_eta, rho_zeta, ierr)
    end subroutine laplace_density

    pure subroutine laplace_solve(ni, nj, nk, x, y, z, phi, m_in, gamma, &
                                   omega, omega_rho, density_update_stride, &
                                   max_iter, tol, &
                                   n_iter_done, ierr, &
                                   progress_iter, max_r_history, rms_r_history, stop_flag)
        integer(ik), intent(in)                                    :: ni, nj, nk
        real(rk),    intent(in),    dimension(nj, nk, ni)          :: x, y, z
        real(rk),    intent(inout), dimension(nj, nk, ni)          :: phi
        real(rk),    intent(in)                                    :: m_in, gamma, omega, omega_rho, tol
        integer(ik), intent(in)                                    :: density_update_stride, max_iter
        integer(ik), intent(out)                                   :: n_iter_done
        integer,     intent(out)                                   :: ierr
        integer(ik), intent(inout), optional                       :: progress_iter
        real(rk),    intent(inout), optional, dimension(max_iter)  :: max_r_history, rms_r_history
        integer(ik), intent(in),    optional                       :: stop_flag
        call laplace_core_solve(ni, nj, nk, x, y, z, phi, m_in, gamma, &
                                 omega, omega_rho, density_update_stride, &
                                 max_iter, tol, &
                                 n_iter_done, ierr, &
                                 progress_iter, max_r_history, rms_r_history, stop_flag)
    end subroutine laplace_solve

    pure subroutine laplace_output(ni, nj, nk, x, y, z, phi, gamma, u, v, w, rho, p, t, mach, ierr)
        integer(ik), intent(in)                         :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk, ni) :: x, y, z, phi
        real(rk),    intent(in)                         :: gamma
        real(rk),    intent(out), dimension(nj, nk, ni) :: u, v, w, rho, p, t, mach
        integer,     intent(out)                        :: ierr
        call laplace_core_output(ni, nj, nk, x, y, z, phi, gamma, u, v, w, rho, p, t, mach, ierr)
    end subroutine laplace_output

end module laplace_f_api
