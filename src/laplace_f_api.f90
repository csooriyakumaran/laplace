module laplace_f_api

    use :: laplace_types, only : rk, ik
    use :: laplace_core,  only : laplace_core_kernel
    use :: laplace_core,  only : laplace_core_grid
    use :: laplace_core,  only : laplace_core_metrics_xi
    use :: laplace_core,  only : laplace_core_metrics_eta
    use :: laplace_core,  only : laplace_core_metrics_zeta
    use :: laplace_core,  only : laplace_core_residual
    use :: laplace_core,  only : laplace_core_density

    implicit none

    private

    public :: laplace_kernel
    public :: laplace_grid
    public :: laplace_metrics_xi
    public :: laplace_metrics_eta
    public :: laplace_metrics_zeta
    public :: laplace_residual
    public :: laplace_density

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

    pure subroutine laplace_residual(ni, nj, nk, phi, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, ierr)
        integer(ik), intent(in)                               :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj,   nk,   ni)   :: phi
        real(rk),    intent(in),  dimension(nj,   nk,   ni-1) :: Aii, Aij, Aik
        real(rk),    intent(in),  dimension(nj-1, nk,   ni)   :: Aji, Ajj, Ajk
        real(rk),    intent(in),  dimension(nj,   nk-1, ni)   :: Aki, Akj, Akk
        real(rk),    intent(out), dimension(nj-2, nk-2, ni-2) :: r
        integer,     intent(out)                              :: ierr
        call laplace_core_residual(ni, nj, nk, phi, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, ierr)
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

end module laplace_f_api
