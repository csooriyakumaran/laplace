module laplace_core
    use :: laplace_types, only : rk, ik

    implicit none

    private
    public :: laplace_core_kernel       ! placeholder to test built system todo(chris): remove this
    public :: laplace_core_grid         ! generates a tranformed grid with tanh clustering near the walls
    public :: laplace_core_metrics_xi   ! face metrics
    public :: laplace_core_metrics_eta  ! face metrics
    public :: laplace_core_metrics_zeta ! face metrics
    public :: laplace_core_residual     ! discrete residual at interior nodes
    public :: laplace_core_density      ! computes density at faces


contains

    !> *******************************************************************************
    !! *  kernel
    !! *
    !! *   Reference: 
    !! *
    !! *   @param p1       Parameter 1
    !! *
    !! *   @return p4      Returns
    !! *
    !! *******************************************************************************/
    subroutine laplace_core_kernel(n, x, y)
        integer(ik), intent(in)                 :: n
        real(rk),    intent(in),  dimension (n) :: x
        real(rk),    intent(out), dimension (n) :: y

        integer(ik) :: i

        do i = 1, n
            write(*,*) i, x(i)
            y(i) = x(i) * 2.0
        end do

    end subroutine laplace_core_kernel

    !> *******************************************************************************
    !! *  grid (PUBLIC)
    !! *
    !! *   Reference: 
    !! *
    !! *   Builds the algebraic grid node coordinates from  the streamwise station 
    !! *   distribution and wall contours. Cross-stream coordinates use tanh clustering
    !! *   from the symmetry plane (s=0) to the wall (s=1)
    !! *
    !! *   @param ni, nj, nk    Grid dimensions (streamwise, vertical, spanwise)
    !! *
    !! *   @param xs            Streamwise station coordinates, strictly increaseing (ni)
    !! *
    !! *   @param hs            Wall half-height at each station (ni)
    !! *
    !! *   @param bs            Wall half-width at each station (ni)
    !! *
    !! *   @param beta_j        tanh clustering strength, vertical direction (> 0)
    !! *
    !! *   @param beta_k        tanh clustering strength, spanwise direction (> 0)
    !! *
    !! *   @param x, y, z       Node coordinates, caller-allocated (nj, nk, ni)
    !! *
    !! *   @param ierr          0 = OK; -1 = invalid beta; -2 = non-monotonic xs
    !! *
   !! *******************************************************************************/
    pure subroutine laplace_core_grid(ni, nj, nk, xs, hs, bs, beta_j, beta_k, x, y, z, ierr)
        integer(ik), intent(in)                         :: ni, nj, nk
        real(rk),    intent(in),  dimension(ni)         :: xs, hs, bs
        real(rk),    intent(in)                         :: beta_j, beta_k
        real(rk),    intent(out), dimension(nj, nk, ni) :: x, y, z
        integer,     intent(out)                        :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: eta, s(nj), t(nk)

        ierr = 0

        if (beta_j <= 0.0_rk .or. beta_k <= 0.0_rk) then
            ierr = -1
            return
        end if

        do i = 2, ni
            if (xs(i) <= xs(i - 1)) then
                ierr = -2
                return
            end if
        end do

        ! cross-stream clustering, each computed once and resused across all (i, k) or (i, j)
        do j = 1, nj
            eta  = real(j - 1, rk) / real(nj - 1, rk)
            s(j) = tanh(beta_j * eta) / tanh(beta_j)
        end do

        do k = 1, nk
            eta  = real(k - 1, rk) / real(nk - 1, rk)
            t(k) = tanh(beta_k * eta) / tanh(beta_k)
        end do

        do i = 1, ni
            do k = 1, nk
                do j = 1, nj
                    x(j, k, i) = xs(i)
                    y(j, k, i) = hs(i) * s(j)
                    z(j, k, i) = bs(i) * t(k)
                end do
            end do
        end do

    end subroutine laplace_core_grid

    pure function central_diff(fm, fp) result(d)
        !> (PRIVATE) Central difference w.r.t. computational index: (fp - fm) / 2.
        !! Adjacent computational indices (xi, eta, zeta) are exactly unit-spaced
        !! by construction, so the divisor is always 2 -- never a coordinate value
        real(rk), intent(in) :: fm, fp
        real(rk)             :: d
        d = (fp - fm) / 2.0_rk
    end function central_diff

    pure function fwd_diff(f0, fp) result(d)
        !> (PRIVATE) Forward difference w.r.t. computational index: fp - f0
        real(rk), intent(in) :: f0, fp
        real(rk)             :: d
        d = fp - f0
    end function fwd_diff

    pure function bwd_diff(fm, f0) result(d)
        !> (PRIVATE) Backward difference w.r.t. computational index: f0 - fm
        real(rk), intent(in) :: fm, f0
        real(rk)             :: d
        d = f0 - fm
    end function bwd_diff

    !> *************************************************************************
    !! * xi_face_deriv (PRIVATE)
    !! *
    !! *   Computes the three face derivatives (d_xi, d_eta, d_zeta) of field f 
    !! *   at the xi-face (i+1/2, j, k). The through-direction (xi) is a plain 
    !! *   difference across the face; each cross-direction (eta, zeta) is a 
    !! *   central difference averaged over the two xi-planes bracketing the face,
    !! *   falling back to a one-sided difference at a j/k boundary.
    !! *
    !! *   Preconditions (not checked here -- this is a private senticle primitive,
    !! *   trusted by its caller, not an ABI boundary): i must satisty 1 <= i < ni
    !! *   so that i+1 is a valid plane; nj >= 2 and nk >= 2 so a cross-direction 
    !! *   neighbour always exists on at least one side
    !! *
    !! *   @param f             field being differenced (e.g, x,y,z, or phi)
    !! *
    !! *   @param ni, nj, nk    field dimensions
    !! *
    !! *   @param i, j, k       face location (i+1/2, j, k)
    !! *
    !! *   @param d_xi, d_eta, d_zeta   te three face derivatives
    !! *
    !! ************************************************************************/
    pure subroutine xi_face_deriv(f, ni, nj, nk, i, j, k, d_xi, d_eta, d_zeta)
        integer(ik), intent(in)                        :: ni, nj, nk
        real(rk),    intent(in), dimension(nj, nk, ni) :: f
        integer(ik), intent(in)                        :: i, j, k
        real(rk),    intent(out)                       :: d_xi, d_eta, d_zeta

        ! through direction: face lies exactly between i and i+1, always available
        d_xi = fwd_diff( f(j,k,i), f(j, k, i+1) )

        ! cross direction (eta), averaged over the two bracketing xi_planes
        if (j > 1 .and. j < nj) then
            d_eta = 0.5_rk * ( central_diff( f(j-1, k, i+1), f(j+1, k, i+1) ) &
                             + central_diff( f(j-1, k, i),   f(j+1, k, i) ) )
        else if (j < nj) then
            d_eta = 0.5_rk * ( fwd_diff( f(j, k, i+1), f(j+1, k, i+1)) &
                             + fwd_diff( f(j, k, i),   f(j+1, k, i)))
        else
            d_eta = 0.5_rk * ( bwd_diff( f(j-1, k, i+1), f(j, k, i+1) ) &
                             + bwd_diff( f(j-1, k, i),   f(j, k, i) ) )
        end if

        ! cross direction (zeta), averaged over the two bracketing xi-planes
        if (k > 1 .and. k < nk) then
            d_zeta = 0.5_rk * ( central_diff( f(j, k-1, i+1), f(j, k+1, i+1) ) &
                              + central_diff( f(j, k-1, i),   f(j, k+1, i) ) )
        else if (k < nk) then
            d_zeta = 0.5_rk * ( fwd_diff( f(j, k, i+1), f(j, k+1, i+1) ) &
                              + fwd_diff( f(j, k, i),   f(j, k+1, i) ) )
        else
            d_zeta = 0.5_rk * ( bwd_diff( f(j, k-1, i+1), f(j, k, i+1) ) &
                              + bwd_diff( f(j, k-1, i),   f(j, k, i) ) )
        end if

    end subroutine xi_face_deriv

    !> *************************************************************************
    !! * eta_face_deriv (PRIVATE)
    !! *
    !! *   @param f             field being differenced (e.g, x,y,z, or phi)
    !! *
    !! *   @param ni, nj, nk    field dimensions
    !! *
    !! *   @param i, j, k       face location (i, j+1/2, k)
    !! *
    !! *   @param d_xi, d_eta, d_zeta   te three face derivatives
    !! *
    !! ************************************************************************/
    pure subroutine eta_face_deriv(f, ni, nj, nk, i, j, k, d_xi, d_eta, d_zeta)
        integer(ik), intent(in)                        :: ni, nj, nk
        real(rk),    intent(in), dimension(nj, nk, ni) :: f
        integer(ik), intent(in)                        :: i, j, k
        real(rk),    intent(out)                       :: d_xi, d_eta, d_zeta

        ! through direction: face lies exactly between j and j+1, always available
        d_eta = fwd_diff(f(j, k, i), f(j+1, k, i))

        ! cross direction (xi), averaged over the two bracketing eta-planes
        if (i > 1 .and. i < ni) then
            d_xi = 0.5_rk * ( central_diff(f(j+1, k, i-1), f(j+1, k, i+1) ) &
                            + central_diff(f(j,   k, i-1), f(j,   k, i+1) ) )
        else if (i < ni) then
            d_xi = 0.5_rk * ( fwd_diff( f(j+1, k, i), f(j+1, k, i+1) ) &
                            + fwd_diff( f(j,   k, i), f(j,   k, i+1) ) )
        else
            d_xi = 0.5_rk * ( bwd_diff( f(j+1, k, i-1), f(j+1, k, i) ) &
                            + bwd_diff( f(j,   k, i-1), f(j,   k, i) ) )
        end if

        ! cross direction (zeta), averaged over the two bracketing eta-planes
        if (k > 1 .and. k < nk) then
            d_zeta = 0.5_rk * ( central_diff( f(j+1, k-1, i), f(j+1, k+1, i) ) &
                              + central_diff( f(j,   k-1, i), f(j,   k+1, i) ) )
        else if (k < nk) then
            d_zeta = 0.5_rk * ( fwd_diff( f(j+1, k, i), f(j+1, k+1, i) ) &
                              + fwd_diff( f(j,   k, i), f(j,   k+1, i) ) )
        else
            d_zeta = 0.5_rk * ( bwd_diff( f(j+1, k-1, i), f(j+1, k, i) ) &
                              + bwd_diff( f(j,   k-1, i), f(j,   k, i) ) )
        end if

    end subroutine eta_face_deriv

    !> *************************************************************************
    !! * zeta_face_deriv (PRIVATE)
    !! *
    !! *   @param f             field being differenced (e.g, x,y,z, or phi)
    !! *
    !! *   @param ni, nj, nk    field dimensions
    !! *
    !! *   @param i, j, k       face location (i, j, k+1/2)
    !! *
    !! *   @param d_xi, d_eta, d_zeta   te three face derivatives
    !! *
    !! ************************************************************************/
    pure subroutine zeta_face_deriv(f, ni, nj, nk, i, j, k, d_xi, d_eta, d_zeta)
        integer(ik), intent(in)                        :: ni, nj, nk
        real(rk),    intent(in), dimension(nj, nk, ni) :: f
        integer(ik), intent(in)                        :: i, j, k
        real(rk),    intent(out)                       :: d_xi, d_eta, d_zeta

        ! through direction: face lies exactly between k and k+1, always available
        d_zeta = fwd_diff(f(j, k, i), f(j, k+1, i))

        ! cross direction (xi), averaged over the two bracketing zeta-planes
        if (i > 1 .and. i < ni) then
            d_xi = 0.5_rk * ( central_diff(f(j, k+1, i-1), f(j, k+1, i+1) ) &
                            + central_diff(f(j, k,   i-1), f(j, k,   i+1) ) )
        else if (i < ni) then
            d_xi = 0.5_rk * ( fwd_diff( f(j, k+1, i), f(j, k+1, i+1) ) &
                            + fwd_diff( f(j, k,   i), f(j, k,   i+1) ) )
        else
            d_xi = 0.5_rk * ( bwd_diff( f(j, k+1, i-1), f(j, k+1, i) ) &
                            + bwd_diff( f(j, k,   i-1), f(j, k,   i) ) )
        end if

        ! cross direction (eta), averaged over the two bracketing zeta-planes
        if (j > 1 .and. j < nj) then
            d_eta = 0.5_rk * ( central_diff( f(j-1, k+1, i), f(j+1, k+1, i) ) &
                             + central_diff( f(j-1, k,   i), f(j+1, k,   i) ) )
        else if (j < nj) then
            d_eta = 0.5_rk * ( fwd_diff( f(j, k+1, i), f(j+1, k+1, i) ) &
                             + fwd_diff( f(j, k,   i), f(j+1, k,   i) ) )
        else
            d_eta = 0.5_rk * ( bwd_diff( f(j-1, k+1, i), f(j, k+1, i) ) &
                             + bwd_diff( f(j-1, k,   i), f(j, k,   i) ) )
        end if

    end subroutine zeta_face_deriv

    !> *************************************************************************
    !! * laplace_core_metrics_xi (PUBLIC)
    !! *
    !! *   Computes the three xi-face metric coefficients from the grid node
    !! *   coordinates and the current face density, evaluated at (i+1/2, j, k).
    !! *
    !! *   @param ni, nj, nk    Grid dimesnions (streamwise, vertical, spanwise)
    !! *
    !! *   @param x, y, z       Node coordinates (nj, nk, ni), e.g. from
    !! *                        laplace_core_grid
    !! *
    !! *   @param rho_xi        Nondimensional density, rho/rho_0 at each xi-face;
    !! *                        caller-supplied, updated during the solve (nj, nk,
    !! *                        ni-1)
    !! *
    !! *   @param Aii, Aij, Aik    A^xixi, A^xieta, A^xizeta at each
    !! *                           xi-face (nj, nk, ni-1)
    !! *
    !! *   @param ierr          0 = OK; -1 - degenerate grid (x_xi <= 0)
    !! *
    !! ************************************************************************/
    pure subroutine laplace_core_metrics_xi(ni, nj, nk, x, y, z, rho_xi, Aii, Aij, Aik, ierr)
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk, ni)   :: x, y, z
        real(rk),    intent(in),  dimension(nj, nk, ni-1) :: rho_xi
        real(rk),    intent(out), dimension(nj, nk, ni-1) :: Aii, Aij, Aik
        integer,     intent(out)                          :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: x_xi, x_eta, x_zeta
        real(rk)    :: y_xi, y_eta, y_zeta
        real(rk)    :: z_xi, z_eta, z_zeta

        ierr = 0

        do i = 1, ni - 1
            do k = 1, nk
                do j = 1, nj
                    call xi_face_deriv(x, ni, nj, nk, i, j, k, x_xi, x_eta, x_zeta)
                    call xi_face_deriv(y, ni, nj, nk, i, j, k, y_xi, y_eta, y_zeta)
                    call xi_face_deriv(z, ni, nj, nk, i, j, k, z_xi, z_eta, z_zeta)

                    if (x_xi <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    Aii(j, k, i) =  rho_xi(j, k, i) * y_eta * z_zeta / x_xi
                    Aij(j, k, i) = -rho_xi(j, k, i) * y_xi  * z_zeta / x_xi
                    Aik(j, k, i) = -rho_xi(j, k, i) * z_xi  * y_eta  / x_xi
                end do
            end do
        end do

    end subroutine laplace_core_metrics_xi

    !> *************************************************************************
    !! * laplace_core_metrics_eta (PUBLIC)
    !! *
    !! *   Computes the three eta-face metric coefficients from the grid node
    !! *   coordinates and the current face density, evaluated at (i, j+1/2, k).
    !! *
    !! *   @param ni, nj, nk    Grid dimesnions (streamwise, vertical, spanwise)
    !! *
    !! *   @param x, y, z       Node coordinates (nj, nk, ni), e.g. from 
    !! *                        laplace_core_grid
    !! *
    !! *   @param rho_eta       Nondimensional density, rho/rho_0 at each eta-face;
    !! *                        caller-supplied, updated during the solve (nj-1,
    !! *                        nk, ni)
    !! *
    !! *   @param Aji, Ajj, Ajk    A^etaxi, A^etaeta, A^etazeta at
    !! *                           each eta-face (nj-1, nk, ni)
    !! *
    !! *   @param ierr          0 = OK; -1 - degenerate grid (x_xi <= 0 or y_eta <= 0)
    !! *
    !! ************************************************************************/
    pure subroutine laplace_core_metrics_eta(ni, nj, nk, x, y, z, rho_eta, Aji, Ajj, Ajk, ierr)
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj,   nk, ni) :: x, y, z
        real(rk),    intent(in),  dimension(nj-1, nk, ni) :: rho_eta
        real(rk),    intent(out), dimension(nj-1, nk, ni) :: Aji, Ajj, Ajk
        integer,     intent(out)                          :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: x_xi, x_eta, x_zeta
        real(rk)    :: y_xi, y_eta, y_zeta
        real(rk)    :: z_xi, z_eta, z_zeta

        ierr = 0

        do i = 1, ni
            do k = 1, nk
                do j = 1, nj -1
                    call eta_face_deriv(x, ni, nj, nk, i, j, k, x_xi, x_eta, x_zeta)
                    call eta_face_deriv(y, ni, nj, nk, i, j, k, y_xi, y_eta, y_zeta)
                    call eta_face_deriv(z, ni, nj, nk, i, j, k, z_xi, z_eta, z_zeta)

                    if (x_xi <= 0.0_rk .or. y_eta <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    Aji(j, k, i) = -rho_eta(j, k, i) * y_xi * z_zeta / x_xi
                    Ajj(j, k, i) =  rho_eta(j, k, i) * x_xi * z_zeta / y_eta * (1.0_rk + (y_xi / x_xi)**2)
                    Ajk(j, k, i) =  rho_eta(j, k, i) * y_xi * z_xi   / x_xi
                end do
            end do
        end do

    end subroutine laplace_core_metrics_eta

    !> *************************************************************************
    !! * laplace_core_metrics_zeta (PUBLIC)
    !! *
    !! *   Computes the three zeta-face metric coefficients from the grid node
    !! *   coordinates and the current face density, evaluated at (i, j, k+1/2).
    !! *
    !! *   @param ni, nj, nk    Grid dimesnions (streamwise, vertical, spanwise)
    !! *
    !! *   @param x, y, z       Node coordinates (nj, nk, ni), e.g. from 
    !! *                        laplace_core_grid
    !! *
    !! *   @param rho_zeta      Nondimensional density, rho/rho_0 at each zeta-face;
    !! *                        caller-supplied, updated during the solve (nj,
    !! *                        nk-1, ni)
    !! *
    !! *   @param Aki, Akj, Akk    A^zetaxi, A^zetaeta, A^zetazeta
    !! *                           at each zeta-face (nj, nk-1, ni)
    !! *
    !! *   @param ierr          0 = OK; -1 - degenerate grid (x_xi <= 0 or z_zeta <= 0)
    !! *
    !! ************************************************************************/
    pure subroutine laplace_core_metrics_zeta(ni, nj, nk, x, y, z, rho_zeta, Aki, Akj, Akk, ierr)
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk,   ni) :: x, y, z
        real(rk),    intent(in),  dimension(nj, nk-1, ni) :: rho_zeta
        real(rk),    intent(out), dimension(nj, nk-1, ni) :: Aki, Akj, Akk
        integer,     intent(out)                          :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: x_xi, x_eta, x_zeta
        real(rk)    :: y_xi, y_eta, y_zeta
        real(rk)    :: z_xi, z_eta, z_zeta

        ierr = 0

        do i = 1, ni
            do k = 1, nk - 1
                do j = 1, nj
                    call zeta_face_deriv(x, ni, nj, nk, i, j, k, x_xi, x_eta, x_zeta)
                    call zeta_face_deriv(y, ni, nj, nk, i, j, k, y_xi, y_eta, y_zeta)
                    call zeta_face_deriv(z, ni, nj, nk, i, j, k, z_xi, z_eta, z_zeta)

                    if (x_xi <= 0.0_rk .or. z_zeta <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    Aki(j, k, i) = -rho_zeta(j, k, i) * z_xi * y_eta / x_xi
                    Akj(j, k, i) =  rho_zeta(j, k, i) * y_xi * z_xi  / x_xi
                    Akk(j, k, i) =  rho_zeta(j, k, i) * x_xi * y_eta / z_zeta * (1.0_rk + (z_xi / x_xi)**2)
                end do
            end do
        end do

    end subroutine laplace_core_metrics_zeta

    !> *************************************************************************
    !! * laplace_core_residual (PUBLIC)
    !! *
    !! *   Computes the discrete residual (README §5) at every interior node
    !! *   from phi and the nine face metric coefficients. Boundary conditions
    !! *   (README §6) are not yet implemented, so this only covers nodes whose
    !! *   six bracketing faces are all ordinary interior faces: i = 2..ni-1,
    !! *   j = 2..nj-1, k = 2..nk-1.
    !! *
    !! *   @param ni, nj, nk              Grid dimensions
    !! *
    !! *   @param phi                     Potential field (nj, nk, ni)
    !! *
    !! *   @param Aii, Aij, Aik           Metrics at xi-faces (nj, nk, ni-1),
    !! *                                  from laplace_core_metrics_xi
    !! *
    !! *   @param Aji, Ajj, Ajk           Metrics at eta-faces (nj-1, nk, ni),
    !! *                                  from laplace_core_metrics_eta
    !! *
    !! *   @param Aki, Akj, Akk           Metrics at zeta-faces (nj, nk-1, ni),
    !! *                                  from laplace_core_metrics_zeta
    !! *
    !! *   @param r                       Residual at each interior node
    !! *                                  (nj-2, nk-2, ni-2); r(jj,kk,ii) is the
    !! *                                  residual at node (ii+1, jj+1, kk+1)
    !! *
    !! *   @param ierr                    0 = OK; -1 = grid too small to have
    !! *                                  any interior node (ni, nj or nk < 3)
    !! *
    !! ************************************************************************/
    pure subroutine laplace_core_residual(ni, nj, nk, phi, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, ierr)
        integer(ik), intent(in)                               :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj,   nk,   ni)   :: phi
        real(rk),    intent(in),  dimension(nj,   nk,   ni-1) :: Aii, Aij, Aik
        real(rk),    intent(in),  dimension(nj-1, nk,   ni)   :: Aji, Ajj, Ajk
        real(rk),    intent(in),  dimension(nj,   nk-1, ni)   :: Aki, Akj, Akk
        real(rk),    intent(out), dimension(nj-2, nk-2, ni-2) :: r
        integer,     intent(out)                              :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: phi_i, phi_j, phi_k
        real(rk)    :: F_p, F_m, G_p, G_m, H_p, H_m

        ierr = 0

        if (ni < 3 .or. nj < 3 .or. nk < 3) then
            ierr = -1
            return
        end if

        do i = 2, ni - 1
            do k = 2, nk - 1
                do j = 2, nj - 1

                    call xi_face_deriv(phi, ni, nj, nk, i, j, k, phi_i, phi_j, phi_k)
                    F_p = Aii(j, k, i) * phi_i + Aij(j, k, i) * phi_j + Aik(j, k, i) * phi_k

                    call xi_face_deriv(phi, ni, nj, nk, i-1, j, k, phi_i, phi_j, phi_k)
                    F_m = Aii(j, k, i-1) * phi_i + Aij(j, k, i-1) * phi_j + Aik(j, k, i-1) * phi_k

                    call eta_face_deriv(phi, ni, nj, nk, i, j, k, phi_i, phi_j, phi_k)
                    G_p = Aji(j, k, i) * phi_i + Ajj(j, k, i) * phi_j + Ajk(j, k, i) * phi_k

                    call eta_face_deriv(phi, ni, nj, nk, i, j-1, k, phi_i, phi_j, phi_k)
                    G_m = Aji(j-1, k, i) * phi_i + Ajj(j-1, k, i) * phi_j + Ajk(j-1, k, i) * phi_k

                    call zeta_face_deriv(phi, ni, nj, nk, i, j, k, phi_i, phi_j, phi_k)
                    H_p = Aki(j, k, i) * phi_i + Akj(j, k, i) * phi_j + Akk(j, k, i) * phi_k

                    call zeta_face_deriv(phi, ni, nj, nk, i, j, k-1, phi_i, phi_j, phi_k)
                    H_m = Aki(j, k-1, i) * phi_i + Akj(j, k-1, i) * phi_j + Akk(j, k-1, i) * phi_k

                    r(j-1, k-1, i-1) = (F_p - F_m) + (G_p - G_m) + (H_p - H_m)

                end do
            end do
        end do

    end subroutine laplace_core_residual

    pure subroutine physical_velocity(x_xi, y_xi, y_eta, z_xi, z_zeta, phi_xi, phi_eta, phi_zeta, u, v, w)
        !> (PRIVATE) Physical velocity from phi's derivatives and the raw metric
        !! terms. Works at any location -- a face (for density closure), or node (for output).
        real(rk), intent(in)  :: x_xi, y_xi, y_eta, z_xi, z_zeta
        real(rk), intent(in)  :: phi_xi, phi_eta, phi_zeta
        real(rk), intent(out) :: u, v, w

        u = phi_xi   / x_xi - (y_xi / (x_xi * y_eta)) * phi_eta - (z_xi / (x_xi * z_zeta)) * phi_zeta
        v = phi_eta  / y_eta
        w = phi_zeta / z_zeta
    end subroutine physical_velocity

    !> *************************************************************************
    !! * laplace_core_density (PUBLIC)
    !! *
    !! *   Computes the nondimensional density rho/rho_0 (README Sec1, Sec2) at
    !! *   each xi-, eta-, and zeta-face from the current phi and the grid, via
    !! *   the physical velocity (README Sec4) and isentropic closure. Mirrors
    !! *   metrics_xi/eta/zeta's own face-gathering, re-deriving the same raw
    !! *   metric terms rather than threading them through from metrics.
    !! *
    !! *   @param ni, nj, nk              Grid dimensions (streamwise, vertical, spanwise)
    !! *
    !! *   @param x, y, z                 Node coordinates (nj, nk, ni)
    !! *
    !! *   @param phi                     Potential field (nj, nk, ni)
    !! *
    !! *   @param gamma                   Ratio of specific heats (> 1)
    !! *
    !! *   @param rho_xi                  Density at each xi-face (nj, nk, ni-1)
    !! *
    !! *   @param rho_eta                 Density at each eta-face (nj-1, nk, ni)
    !! *
    !! *   @param rho_zeta                Density at each zeta-face (nj, nk-1, ni)
    !! *
    !! *   @param ierr                    0 = OK; -1 = degenerate grid
    !! *                                  (x_xi, y_eta or z_zeta <= 0); -2 =
    !! *                                  sonic limit reached (q^2 >= 2/(gamma+1));
    !! *                                  -3 = invalid gamma (<= 1)
    !! *
    !! ************************************************************************/
    pure subroutine laplace_core_density(ni, nj, nk, x, y, z, phi, gamma, rho_xi, rho_eta, rho_zeta, ierr)
        integer(ik), intent(in)                             :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk, ni)     :: x, y, z, phi
        real(rk),    intent(in)                             :: gamma
        real(rk),    intent(out), dimension(nj,   nk, ni-1) :: rho_xi
        real(rk),    intent(out), dimension(nj-1, nk,   ni) :: rho_eta
        real(rk),    intent(out), dimension(nj,   nk-1, ni) :: rho_zeta
        integer,     intent(out)                            :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: x_xi, x_eta, x_zeta
        real(rk)    :: y_xi, y_eta, y_zeta
        real(rk)    :: z_xi, z_eta, z_zeta
        real(rk)    :: phi_xi, phi_eta, phi_zeta
        real(rk)    :: u, v, w, q2, a2, q2_sonic

        ierr = 0

        if (gamma <= 1.0_rk) then
            ierr = -3
            return
        end if

        q2_sonic = 2.0_rk / (gamma + 1.0_rk)

        ! xi-faces
        do i = 1, ni - 1
            do k = 1, nk
                do j = 1, nj
                    call xi_face_deriv(x,   ni, nj, nk, i, j, k, x_xi, x_eta, x_zeta)
                    call xi_face_deriv(y,   ni, nj, nk, i, j, k, y_xi, y_eta, y_zeta)
                    call xi_face_deriv(z,   ni, nj, nk, i, j, k, z_xi, z_eta, z_zeta)
                    call xi_face_deriv(phi, ni, nj, nk, i, j, k, phi_xi, phi_eta, phi_zeta)

                    if (x_xi <= 0.0_rk .or. y_eta <= 0.0_rk .or. z_zeta <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    call physical_velocity(x_xi, y_xi, y_eta, z_xi, z_zeta, phi_xi, phi_eta, phi_zeta, u, v, w)

                    q2 = u*u + v*v + w*w
                    if (q2 >= q2_sonic) then
                        ierr = -2
                        return
                    end if

                    a2 = 1.0_rk - 0.5_rk * (gamma - 1.0_rk) * q2
                    rho_xi(j, k, i) = a2 ** (1.0_rk / (gamma - 1.0_rk))
                end do
            end do
        end do

        ! eta-faces
        do i = 1, ni
            do k = 1, nk
                do j = 1, nj - 1
                    call eta_face_deriv(x,   ni, nj, nk, i, j, k, x_xi, x_eta, x_zeta)
                    call eta_face_deriv(y,   ni, nj, nk, i, j, k, y_xi, y_eta, y_zeta)
                    call eta_face_deriv(z,   ni, nj, nk, i, j, k, z_xi, z_eta, z_zeta)
                    call eta_face_deriv(phi, ni, nj, nk, i, j, k, phi_xi, phi_eta, phi_zeta)

                    if (x_xi <= 0.0_rk .or. y_eta <= 0.0_rk .or. z_zeta <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    call physical_velocity(x_xi, y_xi, y_eta, z_xi, z_zeta, phi_xi, phi_eta, phi_zeta, u, v, w)

                    q2 = u*u + v*v + w*w
                    if (q2 >= q2_sonic) then
                        ierr = -2
                        return
                    end if

                    a2 = 1.0_rk - 0.5_rk * (gamma - 1.0_rk) * q2
                    rho_eta(j, k, i) = a2 ** (1.0_rk / (gamma - 1.0_rk))
                end do
            end do
        end do

        ! zeta-faces
        do i = 1, ni
            do k = 1, nk - 1
                do j = 1, nj
                    call zeta_face_deriv(x,   ni, nj, nk, i, j, k, x_xi, x_eta, x_zeta)
                    call zeta_face_deriv(y,   ni, nj, nk, i, j, k, y_xi, y_eta, y_zeta)
                    call zeta_face_deriv(z,   ni, nj, nk, i, j, k, z_xi, z_eta, z_zeta)
                    call zeta_face_deriv(phi, ni, nj, nk, i, j, k, phi_xi, phi_eta, phi_zeta)

                    if (x_xi <= 0.0_rk .or. y_eta <= 0.0_rk .or. z_zeta <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    call physical_velocity(x_xi, y_xi, y_eta, z_xi, z_zeta, phi_xi, phi_eta, phi_zeta, u, v, w)

                    q2 = u*u + v*v + w*w
                    if (q2 >= q2_sonic) then
                        ierr = -2
                        return
                    end if

                    a2 = 1.0_rk - 0.5_rk * (gamma - 1.0_rk) * q2
                    rho_zeta(j, k, i) = a2 ** (1.0_rk / (gamma - 1.0_rk))
                end do
            end do
        end do
    end subroutine laplace_core_density

end module laplace_core
