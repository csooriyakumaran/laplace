module laplace_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int64_t, c_int32_t
    use            :: laplace_types, only : rk, ik
    use            :: laplace_core,  only : laplace_core_kernel
    use            :: laplace_core,  only : laplace_core_grid
    use            :: laplace_core,  only : laplace_core_metrics_xi
    use            :: laplace_core,  only : laplace_core_metrics_eta
    use            :: laplace_core,  only : laplace_core_metrics_zeta
    use            :: laplace_core,  only : laplace_core_residual

    ! add use statement for all public api functions

    implicit none

    private

    public :: laplace_kernel_c
    public :: laplace_grid_c
    public :: laplace_metrics_xi_c
    public :: laplace_metrics_eta_c
    public :: laplace_metrics_zeta_c
    public :: laplace_residual_c

contains

    subroutine laplace_kernel_c(n, x, y) bind(C, name="laplace_kernel")
        integer(ik), intent(in),  value        :: n
        real(rk),    intent(in),  dimension(n) :: x
        real(rk),    intent(out), dimension(n) :: y
        call laplace_core_kernel(n, x, y)
    end subroutine laplace_kernel_c

    pure subroutine laplace_grid_c(ni, nj, nk, xs, hs, bs, beta_j, beta_k, x, y, z, ierr) bind(C, name="laplace_grid")
        integer(ik), intent(in),  value                 :: ni, nj, nk
        real(rk),    intent(in),  dimension(ni)         :: xs, hs, bs
        real(rk),    intent(in),  value                 :: beta_j, beta_k
        real(rk),    intent(out), dimension(nj, nk, ni) :: x, y, z
        integer(c_int32_t), intent(out)                 :: ierr
        call laplace_core_grid(ni, nj, nk, xs, hs, bs, beta_j, beta_k, x, y, z, ierr)
    end subroutine laplace_grid_c

    pure subroutine laplace_metrics_xi_c(ni, nj, nk, x, y, z, rho_xi, Aii, Aij, Aik, ierr) bind(C, name="laplace_metrics_xi")
        integer(ik), intent(in),  value                   :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk, ni)   :: x, y, z
        real(rk),    intent(in),  dimension(nj, nk, ni-1) :: rho_xi
        real(rk),    intent(out), dimension(nj, nk, ni-1) :: Aii, Aij, Aik
        integer(c_int32_t), intent(out)                   :: ierr
        call laplace_core_metrics_xi(ni, nj, nk, x, y, z, rho_xi, Aii, Aij, Aik, ierr)
    end subroutine laplace_metrics_xi_c

    pure subroutine laplace_metrics_eta_c(ni, nj, nk, x, y, z, rho_eta, Aji, Ajj, Ajk, ierr) bind(C, name="laplace_metrics_eta")
        integer(ik), intent(in),  value                   :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj,   nk, ni) :: x, y, z
        real(rk),    intent(in),  dimension(nj-1, nk, ni) :: rho_eta
        real(rk),    intent(out), dimension(nj-1, nk, ni) :: Aji, Ajj, Ajk
        integer(c_int32_t), intent(out)                   :: ierr
        call laplace_core_metrics_eta(ni, nj, nk, x, y, z, rho_eta, Aji, Ajj, Ajk, ierr)
    end subroutine laplace_metrics_eta_c

    pure subroutine laplace_metrics_zeta_c(ni, nj, nk, x, y, z, rho_zeta, Aki, Akj, Akk, ierr) bind(C, name="laplace_metrics_zeta")
        integer(ik), intent(in),  value                   :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk,   ni) :: x, y, z
        real(rk),    intent(in),  dimension(nj, nk-1, ni) :: rho_zeta
        real(rk),    intent(out), dimension(nj, nk-1, ni) :: Aki, Akj, Akk
        integer(c_int32_t), intent(out)                   :: ierr
        call laplace_core_metrics_zeta(ni, nj, nk, x, y, z, rho_zeta, Aki, Akj, Akk, ierr)
    end subroutine laplace_metrics_zeta_c

    pure subroutine laplace_residual_c(ni, nj, nk, phi, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, ierr) bind(C, name="laplace_residual")
        integer(ik), intent(in),  value                       :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj,   nk,   ni)   :: phi
        real(rk),    intent(in),  dimension(nj,   nk,   ni-1) :: Aii, Aij, Aik
        real(rk),    intent(in),  dimension(nj-1, nk,   ni)   :: Aji, Ajj, Ajk
        real(rk),    intent(in),  dimension(nj,   nk-1, ni)   :: Aki, Akj, Akk
        real(rk),    intent(out), dimension(nj-2, nk-2, ni-2) :: r
        integer(c_int32_t),  intent(out)                   :: ierr
        call laplace_core_residual(ni, nj, nk, phi, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, ierr)
    end subroutine laplace_residual_c

end module laplace_c_api
