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
    public :: laplace_core_solve        ! SLOR convergence loop with density-update cadence
    public :: laplace_core_output       ! node-centered velocity/rho/p/T/M output
    public :: thomas_solve              ! tridiagonal solve (README Sec7); exposed for direct unit testing only, not wrapped by f_api/c_api
    public :: sweep                     ! one SLOR sweep (README Sec7); exposed for direct unit testing only, not wrapped by f_api/c_api
    public :: laplace_atomic_store_rel_i64 ! release-store (src/laplace_atomic.c); exposed for direct unit testing only, not wrapped by f_api/c_api
    public :: laplace_atomic_load_acq_i64  ! acquire-load  (src/laplace_atomic.c); exposed for direct unit testing only, not wrapped by f_api/c_api

    !> PURE interfaces to the two atomic primitives in src/laplace_atomic.c.
    !! PURE here is a programmer assertion, not something the compiler can
    !! verify across the C/Fortran boundary -- the same basis on which any
    !! PURE procedure's call to an external bind(C) routine is trusted. The
    !! actual purpose of these two calls is a controlled, intentional side
    !! effect (an atomic store/load through caller-owned memory), which is
    !! exactly the class of effect PURE already permits via intent(inout)/
    !! intent(out) dummy arguments -- this simply reaches that same effect
    !! through memory Fortran itself has no portable way to mark atomic.
    interface
        pure subroutine laplace_atomic_store_rel_i64(p, v) bind(C, name="laplace_atomic_store_rel_i64")
            import :: ik
            integer(ik), intent(inout)     :: p
            integer(ik), intent(in), value :: v
        end subroutine laplace_atomic_store_rel_i64

        !> witness: pass a value that changes on every call (e.g. a loop
        !! counter) -- its only purpose is to prevent a PURE-licensed
        !! compiler from treating repeated calls as having "the same
        !! arguments" and caching/hoisting the load. See laplace_atomic.c.
        pure function laplace_atomic_load_acq_i64(p, witness) result(v) &
                bind(C, name="laplace_atomic_load_acq_i64")
            import :: ik
            integer(ik), intent(in)        :: p
            integer(ik), intent(in), value :: witness
            integer(ik)                    :: v
        end function laplace_atomic_load_acq_i64
    end interface


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

    pure subroutine raw_metrics_xi(ni, nj, nk, x, y, z, x_xi, y_xi, y_eta, z_xi, z_zeta, ierr)
        !> (PRIVATE) The 5 raw grid-derivative terms needed by both
        !! laplace_core_metrics_xi and laplace_core_density at every xi-face --
        !! depend only on the grid (x, y, z), never on phi or rho. Factored out
        !! so laplace_core_solve can compute these once per solve instead of on
        !! every density/metrics update (the grid never changes mid-solve).
        !! Degenerate-grid check matches laplace_core_metrics_xi's original
        !! (x_xi <= 0 only) -- density's stricter 3-condition check is done by
        !! its own consumer (density_xi_from_raw), not duplicated here.
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk, ni)   :: x, y, z
        real(rk),    intent(out), dimension(nj, nk, ni-1) :: x_xi, y_xi, y_eta, z_xi, z_zeta
        integer,     intent(out)                          :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: xx, xe, xz, yx, ye, yz, zx, ze, zz

        ierr = 0

        do i = 1, ni - 1
            do k = 1, nk
                do j = 1, nj
                    call xi_face_deriv(x, ni, nj, nk, i, j, k, xx, xe, xz)
                    call xi_face_deriv(y, ni, nj, nk, i, j, k, yx, ye, yz)
                    call xi_face_deriv(z, ni, nj, nk, i, j, k, zx, ze, zz)

                    if (xx <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    x_xi(j, k, i) = xx; y_xi(j, k, i) = yx; y_eta(j, k, i) = ye
                    z_xi(j, k, i) = zx; z_zeta(j, k, i) = zz
                end do
            end do
        end do
    end subroutine raw_metrics_xi

    pure subroutine raw_metrics_eta(ni, nj, nk, x, y, z, x_xi, y_xi, y_eta, z_xi, z_zeta, ierr)
        !> (PRIVATE) Same as raw_metrics_xi, evaluated at eta-faces. Check
        !! matches laplace_core_metrics_eta's original (x_xi <= 0 .or. y_eta <= 0).
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj,   nk, ni) :: x, y, z
        real(rk),    intent(out), dimension(nj-1, nk, ni) :: x_xi, y_xi, y_eta, z_xi, z_zeta
        integer,     intent(out)                          :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: xx, xe, xz, yx, ye, yz, zx, ze, zz

        ierr = 0

        do i = 1, ni
            do k = 1, nk
                do j = 1, nj - 1
                    call eta_face_deriv(x, ni, nj, nk, i, j, k, xx, xe, xz)
                    call eta_face_deriv(y, ni, nj, nk, i, j, k, yx, ye, yz)
                    call eta_face_deriv(z, ni, nj, nk, i, j, k, zx, ze, zz)

                    if (xx <= 0.0_rk .or. ye <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    x_xi(j, k, i) = xx; y_xi(j, k, i) = yx; y_eta(j, k, i) = ye
                    z_xi(j, k, i) = zx; z_zeta(j, k, i) = zz
                end do
            end do
        end do
    end subroutine raw_metrics_eta

    pure subroutine raw_metrics_zeta(ni, nj, nk, x, y, z, x_xi, y_xi, y_eta, z_xi, z_zeta, ierr)
        !> (PRIVATE) Same as raw_metrics_xi, evaluated at zeta-faces. Check
        !! matches laplace_core_metrics_zeta's original (x_xi <= 0 .or. z_zeta <= 0).
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk,   ni) :: x, y, z
        real(rk),    intent(out), dimension(nj, nk-1, ni) :: x_xi, y_xi, y_eta, z_xi, z_zeta
        integer,     intent(out)                          :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: xx, xe, xz, yx, ye, yz, zx, ze, zz

        ierr = 0

        do i = 1, ni
            do k = 1, nk - 1
                do j = 1, nj
                    call zeta_face_deriv(x, ni, nj, nk, i, j, k, xx, xe, xz)
                    call zeta_face_deriv(y, ni, nj, nk, i, j, k, yx, ye, yz)
                    call zeta_face_deriv(z, ni, nj, nk, i, j, k, zx, ze, zz)

                    if (xx <= 0.0_rk .or. zz <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    x_xi(j, k, i) = xx; y_xi(j, k, i) = yx; y_eta(j, k, i) = ye
                    z_xi(j, k, i) = zx; z_zeta(j, k, i) = zz
                end do
            end do
        end do
    end subroutine raw_metrics_zeta

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

        real(rk), allocatable, dimension(:,:,:) :: x_xi, y_xi, y_eta, z_xi, z_zeta

        allocate(x_xi(nj,nk,ni-1), y_xi(nj,nk,ni-1), y_eta(nj,nk,ni-1), z_xi(nj,nk,ni-1), z_zeta(nj,nk,ni-1))

        call raw_metrics_xi(ni, nj, nk, x, y, z, x_xi, y_xi, y_eta, z_xi, z_zeta, ierr)
        if (ierr /= 0) return

        Aii =  rho_xi * y_eta * z_zeta / x_xi
        Aij = -rho_xi * y_xi  * z_zeta / x_xi
        Aik = -rho_xi * z_xi  * y_eta  / x_xi

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

        real(rk), allocatable, dimension(:,:,:) :: x_xi, y_xi, y_eta, z_xi, z_zeta

        allocate(x_xi(nj-1,nk,ni), y_xi(nj-1,nk,ni), y_eta(nj-1,nk,ni), z_xi(nj-1,nk,ni), z_zeta(nj-1,nk,ni))

        call raw_metrics_eta(ni, nj, nk, x, y, z, x_xi, y_xi, y_eta, z_xi, z_zeta, ierr)
        if (ierr /= 0) return

        Aji = -rho_eta * y_xi * z_zeta / x_xi
        Ajj =  rho_eta * x_xi * z_zeta / y_eta * (1.0_rk + (y_xi / x_xi)**2)
        Ajk =  rho_eta * y_xi * z_xi   / x_xi

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

        real(rk), allocatable, dimension(:,:,:) :: x_xi, y_xi, y_eta, z_xi, z_zeta

        allocate(x_xi(nj,nk-1,ni), y_xi(nj,nk-1,ni), y_eta(nj,nk-1,ni), z_xi(nj,nk-1,ni), z_zeta(nj,nk-1,ni))

        call raw_metrics_zeta(ni, nj, nk, x, y, z, x_xi, y_xi, y_eta, z_xi, z_zeta, ierr)
        if (ierr /= 0) return

        Aki = -rho_zeta * z_xi * y_eta / x_xi
        Akj =  rho_zeta * y_xi * z_xi  / x_xi
        Akk =  rho_zeta * x_xi * y_eta / z_zeta * (1.0_rk + (z_xi / x_xi)**2)

    end subroutine laplace_core_metrics_zeta

    pure function node_deriv_xi(f, nj, nk, ni, i, j, k) result(d)
        !> (PRIVATE) Central difference of f w.r.t. xi at a single node,
        !! one-sided at an i-boundary. No cross-direction averaging -- unlike
        !! xi_face_deriv, this is for a quantity needed AT a node itself
        !! (node-centered output, README Sec8), not at a face straddling two
        !! planes. See node_deriv_eta/node_deriv_zeta.
        integer(ik), intent(in) :: nj, nk, ni
        real(rk),    intent(in) :: f(nj, nk, ni)
        integer(ik), intent(in) :: i, j, k
        real(rk) :: d

        if (i > 1 .and. i < ni) then
            d = central_diff(f(j, k, i-1), f(j, k, i+1))
        else if (i < ni) then
            d = fwd_diff(f(j, k, i), f(j, k, i+1))
        else
            d = bwd_diff(f(j, k, i-1), f(j, k, i))
        end if
    end function node_deriv_xi

    pure function node_deriv_eta(f, nj, nk, ni, i, j, k) result(d)
        !> (PRIVATE) Central difference of f w.r.t. eta at a single node,
        !! one-sided at a j-boundary. Used both for the inlet boundary's
        !! prescribed mass flux (laplace_core_residual) and for node-centered
        !! output (README Sec8). See node_deriv_xi.
        integer(ik), intent(in) :: nj, nk, ni
        real(rk),    intent(in) :: f(nj, nk, ni)
        integer(ik), intent(in) :: i, j, k
        real(rk) :: d

        if (j > 1 .and. j < nj) then
            d = central_diff(f(j-1, k, i), f(j+1, k, i))
        else if (j < nj) then
            d = fwd_diff(f(j, k, i), f(j+1, k, i))
        else
            d = bwd_diff(f(j-1, k, i), f(j, k, i))
        end if
    end function node_deriv_eta

    pure function node_deriv_zeta(f, nj, nk, ni, i, j, k) result(d)
        !> (PRIVATE) Central difference of f w.r.t. zeta at a single node,
        !! one-sided at a k-boundary. See node_deriv_xi.
        integer(ik), intent(in) :: nj, nk, ni
        real(rk),    intent(in) :: f(nj, nk, ni)
        integer(ik), intent(in) :: i, j, k
        real(rk) :: d

        if (k > 1 .and. k < nk) then
            d = central_diff(f(j, k-1, i), f(j, k+1, i))
        else if (k < nk) then
            d = fwd_diff(f(j, k, i), f(j, k+1, i))
        else
            d = bwd_diff(f(j, k-1, i), f(j, k, i))
        end if
    end function node_deriv_zeta

    pure subroutine node_residual(ni, nj, nk, phi, y, z, m_in, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, i, j, k, r)
        !> (PRIVATE) Computes the discrete residual (README §5) at a single
        !! solved node (i,j,k), 1 <= i <= ni-1, full j, full k. Factored out of
        !! laplace_core_residual so laplace_core_solve can evaluate it one
        !! line at a time instead of recomputing the whole domain per line.
        !! See laplace_core_residual for the README §6 boundary conditions.
        integer(ik), intent(in)                               :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj,   nk,   ni)   :: phi, y, z
        real(rk),    intent(in)                               :: m_in
        real(rk),    intent(in),  dimension(nj,   nk,   ni-1) :: Aii, Aij, Aik
        real(rk),    intent(in),  dimension(nj-1, nk,   ni)   :: Aji, Ajj, Ajk
        real(rk),    intent(in),  dimension(nj,   nk-1, ni)   :: Aki, Akj, Akk
        integer(ik), intent(in)                               :: i, j, k
        real(rk),    intent(out)                              :: r

        real(rk) :: phi_i, phi_j, phi_k
        real(rk) :: F_p, F_m, G_p, G_m, H_p, H_m

        call xi_face_deriv(phi, ni, nj, nk, i, j, k, phi_i, phi_j, phi_k)
        F_p = Aii(j, k, i) * phi_i + Aij(j, k, i) * phi_j + Aik(j, k, i) * phi_k

        if (i > 1) then
            call xi_face_deriv(phi, ni, nj, nk, i-1, j, k, phi_i, phi_j, phi_k)
            F_m = Aii(j, k, i-1) * phi_i + Aij(j, k, i-1) * phi_j + Aik(j, k, i-1) * phi_k
        else
            F_m = m_in * node_deriv_eta(y, nj, nk, ni, i, j, k) &
                       * node_deriv_zeta(z, nj, nk, ni, i, j, k)
        end if

        if (j < nj) then
            call eta_face_deriv(phi, ni, nj, nk, i, j, k, phi_i, phi_j, phi_k)
            G_p = Aji(j, k, i) * phi_i + Ajj(j, k, i) * phi_j + Ajk(j, k, i) * phi_k
        else
            G_p = 0.0_rk
        end if

        if (j > 1) then
            call eta_face_deriv(phi, ni, nj, nk, i, j-1, k, phi_i, phi_j, phi_k)
            G_m = Aji(j-1, k, i) * phi_i + Ajj(j-1, k, i) * phi_j + Ajk(j-1, k, i) * phi_k
        else
            G_m = 0.0_rk
        end if

        if (k < nk) then
            call zeta_face_deriv(phi, ni, nj, nk, i, j, k, phi_i, phi_j, phi_k)
            H_p = Aki(j, k, i) * phi_i + Akj(j, k, i) * phi_j + Akk(j, k, i) * phi_k
        else
            H_p = 0.0_rk
        end if

        if (k > 1) then
            call zeta_face_deriv(phi, ni, nj, nk, i, j, k-1, phi_i, phi_j, phi_k)
            H_m = Aki(j, k-1, i) * phi_i + Akj(j, k-1, i) * phi_j + Akk(j, k-1, i) * phi_k
        else
            H_m = 0.0_rk
        end if

        r = (F_p - F_m) + (G_p - G_m) + (H_p - H_m)

    end subroutine node_residual

    !> *************************************************************************
    !! * laplace_core_residual (PUBLIC)
    !! *
    !! *   Computes the discrete residual (README §5) at every solved node:
    !! *   i = 1..ni-1, full j, full k. The outlet plane (i = ni) is Dirichlet
    !! *   phi = 0 (README §6) -- fixed data, not part of the solved system, so
    !! *   it has no residual of its own; phi(:,:,ni) is read as ordinary known
    !! *   data by the i = ni-1 nodes' F_p.
    !! *
    !! *   Boundary conditions (README §6) folded in directly:
    !! *     - symmetry (j=1, k=1) / wall (j=nj, k=nk): the missing-neighbour
    !! *       face's flux is omitted entirely (zero), not substituted
    !! *     - inlet (i=1): F_m is the prescribed constant mass flux
    !! *       m_in * y_eta * z_zeta, evaluated at the inlet plane directly
    !! *       (not a face average -- there is no i=0 plane to average with)
    !! *
    !! *   @param ni, nj, nk              Grid dimensions
    !! *
    !! *   @param phi                     Potential field (nj, nk, ni)
    !! *
    !! *   @param y, z                    Node coordinates (nj, nk, ni); needed
    !! *                                  only for the inlet's geometric term
    !! *
    !! *   @param m_in                    Prescribed nondimensional inlet mass
    !! *                                  flux, rho*u/(rho_0*a_0) (README §2)
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
    !! *   @param r                       Residual at each solved node
    !! *                                  (nj, nk, ni-1); r(j,k,i) is the
    !! *                                  residual at node (i,j,k) -- no index
    !! *                                  offset, unlike the interior-only version
    !! *
    !! *   @param ierr                    0 = OK; -1 = ni, nj or nk < 2
    !! *
    !! ************************************************************************/
    pure subroutine laplace_core_residual(ni, nj, nk, phi, y, z, m_in, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, ierr)
        integer(ik), intent(in)                               :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj,   nk,   ni)   :: phi, y, z
        real(rk),    intent(in)                               :: m_in
        real(rk),    intent(in),  dimension(nj,   nk,   ni-1) :: Aii, Aij, Aik
        real(rk),    intent(in),  dimension(nj-1, nk,   ni)   :: Aji, Ajj, Ajk
        real(rk),    intent(in),  dimension(nj,   nk-1, ni)   :: Aki, Akj, Akk
        real(rk),    intent(out), dimension(nj,   nk,   ni-1) :: r
        integer,     intent(out)                              :: ierr

        integer(ik) :: i, j, k

        ierr = 0

        if (ni < 2 .or. nj < 2 .or. nk < 2) then
            ierr = -1
            return
        end if

        do i = 1, ni - 1
            do k = 1, nk
                do j = 1, nj
                    call node_residual(ni, nj, nk, phi, y, z, m_in, &
                                        Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, &
                                        i, j, k, r(j, k, i))
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

    pure subroutine density_from_raw(x_xi, y_xi, y_eta, z_xi, z_zeta, phi_xi, phi_eta, phi_zeta, gamma, rho, ierr)
        !> (PRIVATE) Isentropic closure (README Sec1) at a single face, given
        !! already-derived raw metric terms and phi's derivatives there. Shared
        !! by laplace_core_density (fresh raw terms every call) and
        !! laplace_core_solve's hot loop (raw terms cached once per solve,
        !! since the grid never changes).
        real(rk), intent(in)  :: x_xi, y_xi, y_eta, z_xi, z_zeta
        real(rk), intent(in)  :: phi_xi, phi_eta, phi_zeta
        real(rk), intent(in)  :: gamma
        real(rk), intent(out) :: rho
        integer,  intent(out) :: ierr

        real(rk) :: u, v, w, q2, a2, q2_sonic

        ierr = 0

        call physical_velocity(x_xi, y_xi, y_eta, z_xi, z_zeta, phi_xi, phi_eta, phi_zeta, u, v, w)

        q2_sonic = 2.0_rk / (gamma + 1.0_rk)
        q2 = u*u + v*v + w*w
        if (q2 >= q2_sonic) then
            ierr = -2
            return
        end if

        a2  = 1.0_rk - 0.5_rk * (gamma - 1.0_rk) * q2
        rho = a2 ** (1.0_rk / (gamma - 1.0_rk))
    end subroutine density_from_raw

    pure subroutine density_xi_from_raw(ni, nj, nk, phi, x_xi, y_xi, y_eta, z_xi, z_zeta, gamma, rho_xi, ierr)
        !> (PRIVATE) density_from_raw looped over every xi-face, using
        !! caller-supplied raw terms (fresh or cached) instead of re-deriving
        !! x_xi/y_xi/y_eta/z_xi/z_zeta from x,y,z. phi's own derivative is
        !! still derived fresh every call -- phi changes every iteration,
        !! unlike the grid. Degenerate check here matches laplace_core_density's
        !! original (all 3 conditions), stricter than raw_metrics_xi's own.
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk, ni)   :: phi
        real(rk),    intent(in),  dimension(nj, nk, ni-1) :: x_xi, y_xi, y_eta, z_xi, z_zeta
        real(rk),    intent(in)                           :: gamma
        real(rk),    intent(out), dimension(nj, nk, ni-1) :: rho_xi
        integer,     intent(out)                          :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: phi_xi, phi_eta, phi_zeta

        ierr = 0

        do i = 1, ni - 1
            do k = 1, nk
                do j = 1, nj
                    if (x_xi(j,k,i) <= 0.0_rk .or. y_eta(j,k,i) <= 0.0_rk .or. z_zeta(j,k,i) <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    call xi_face_deriv(phi, ni, nj, nk, i, j, k, phi_xi, phi_eta, phi_zeta)

                    call density_from_raw(x_xi(j,k,i), y_xi(j,k,i), y_eta(j,k,i), z_xi(j,k,i), z_zeta(j,k,i), &
                                           phi_xi, phi_eta, phi_zeta, gamma, rho_xi(j,k,i), ierr)
                    if (ierr /= 0) return
                end do
            end do
        end do
    end subroutine density_xi_from_raw

    pure subroutine density_eta_from_raw(ni, nj, nk, phi, x_xi, y_xi, y_eta, z_xi, z_zeta, gamma, rho_eta, ierr)
        !> (PRIVATE) See density_xi_from_raw; evaluated at eta-faces.
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj,   nk, ni) :: phi
        real(rk),    intent(in),  dimension(nj-1, nk, ni) :: x_xi, y_xi, y_eta, z_xi, z_zeta
        real(rk),    intent(in)                           :: gamma
        real(rk),    intent(out), dimension(nj-1, nk, ni) :: rho_eta
        integer,     intent(out)                          :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: phi_xi, phi_eta, phi_zeta

        ierr = 0

        do i = 1, ni
            do k = 1, nk
                do j = 1, nj - 1
                    if (x_xi(j,k,i) <= 0.0_rk .or. y_eta(j,k,i) <= 0.0_rk .or. z_zeta(j,k,i) <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    call eta_face_deriv(phi, ni, nj, nk, i, j, k, phi_xi, phi_eta, phi_zeta)

                    call density_from_raw(x_xi(j,k,i), y_xi(j,k,i), y_eta(j,k,i), z_xi(j,k,i), z_zeta(j,k,i), &
                                           phi_xi, phi_eta, phi_zeta, gamma, rho_eta(j,k,i), ierr)
                    if (ierr /= 0) return
                end do
            end do
        end do
    end subroutine density_eta_from_raw

    pure subroutine density_zeta_from_raw(ni, nj, nk, phi, x_xi, y_xi, y_eta, z_xi, z_zeta, gamma, rho_zeta, ierr)
        !> (PRIVATE) See density_xi_from_raw; evaluated at zeta-faces.
        integer(ik), intent(in)                           :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk,   ni) :: phi
        real(rk),    intent(in),  dimension(nj, nk-1, ni) :: x_xi, y_xi, y_eta, z_xi, z_zeta
        real(rk),    intent(in)                           :: gamma
        real(rk),    intent(out), dimension(nj, nk-1, ni) :: rho_zeta
        integer,     intent(out)                          :: ierr

        integer(ik) :: i, j, k
        real(rk)    :: phi_xi, phi_eta, phi_zeta

        ierr = 0

        do i = 1, ni
            do k = 1, nk - 1
                do j = 1, nj
                    if (x_xi(j,k,i) <= 0.0_rk .or. y_eta(j,k,i) <= 0.0_rk .or. z_zeta(j,k,i) <= 0.0_rk) then
                        ierr = -1
                        return
                    end if

                    call zeta_face_deriv(phi, ni, nj, nk, i, j, k, phi_xi, phi_eta, phi_zeta)

                    call density_from_raw(x_xi(j,k,i), y_xi(j,k,i), y_eta(j,k,i), z_xi(j,k,i), z_zeta(j,k,i), &
                                           phi_xi, phi_eta, phi_zeta, gamma, rho_zeta(j,k,i), ierr)
                    if (ierr /= 0) return
                end do
            end do
        end do
    end subroutine density_zeta_from_raw

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

        real(rk), allocatable, dimension(:,:,:) :: xi_x_xi, xi_y_xi, xi_y_eta, xi_z_xi, xi_z_zeta
        real(rk), allocatable, dimension(:,:,:) :: eta_x_xi, eta_y_xi, eta_y_eta, eta_z_xi, eta_z_zeta
        real(rk), allocatable, dimension(:,:,:) :: zeta_x_xi, zeta_y_xi, zeta_y_eta, zeta_z_xi, zeta_z_zeta

        ierr = 0

        if (gamma <= 1.0_rk) then
            ierr = -3
            return
        end if

        allocate(xi_x_xi(nj,nk,ni-1), xi_y_xi(nj,nk,ni-1), xi_y_eta(nj,nk,ni-1), xi_z_xi(nj,nk,ni-1), xi_z_zeta(nj,nk,ni-1))
        call raw_metrics_xi(ni, nj, nk, x, y, z, xi_x_xi, xi_y_xi, xi_y_eta, xi_z_xi, xi_z_zeta, ierr)
        if (ierr == 0) call density_xi_from_raw(ni, nj, nk, phi, xi_x_xi, xi_y_xi, xi_y_eta, xi_z_xi, xi_z_zeta, gamma, rho_xi, ierr)
        if (ierr /= 0) return

        allocate(eta_x_xi(nj-1,nk,ni), eta_y_xi(nj-1,nk,ni), eta_y_eta(nj-1,nk,ni), eta_z_xi(nj-1,nk,ni), eta_z_zeta(nj-1,nk,ni))
        call raw_metrics_eta(ni, nj, nk, x, y, z, eta_x_xi, eta_y_xi, eta_y_eta, eta_z_xi, eta_z_zeta, ierr)
        if (ierr == 0) call density_eta_from_raw(ni, nj, nk, phi, eta_x_xi, eta_y_xi, eta_y_eta, eta_z_xi, eta_z_zeta, gamma, rho_eta, ierr)
        if (ierr /= 0) return

        allocate(zeta_x_xi(nj,nk-1,ni), zeta_y_xi(nj,nk-1,ni), zeta_y_eta(nj,nk-1,ni), zeta_z_xi(nj,nk-1,ni), zeta_z_zeta(nj,nk-1,ni))
        call raw_metrics_zeta(ni, nj, nk, x, y, z, zeta_x_xi, zeta_y_xi, zeta_y_eta, zeta_z_xi, zeta_z_zeta, ierr)
        if (ierr == 0) call density_zeta_from_raw(ni, nj, nk, phi, zeta_x_xi, zeta_y_xi, zeta_y_eta, zeta_z_xi, zeta_z_zeta, gamma, rho_zeta, ierr)
    end subroutine laplace_core_density

    pure subroutine node_output(ni, nj, nk, x, y, z, phi, gamma, i, j, k, u, v, w, rho, p, t, mach, ierr)
        !> (PRIVATE) Node-centered physical output at a single node (README
        !! Sec8): velocity from node-centered (not face-centered) derivatives
        !! of phi through physical_velocity, then rho/p/T/M from the same
        !! isentropic closure as laplace_core_density. All nondimensional
        !! (README Sec1/Sec2); dimensionalizing by rho_0/p_0/T_0/a_0 is the
        !! caller's job.
        integer(ik), intent(in)                        :: ni, nj, nk
        real(rk),    intent(in), dimension(nj, nk, ni) :: x, y, z, phi
        real(rk),    intent(in)                        :: gamma
        integer(ik), intent(in)                        :: i, j, k
        real(rk),    intent(out)                       :: u, v, w, rho, p, t, mach
        integer,     intent(out)                       :: ierr

        real(rk) :: x_xi, y_xi, y_eta, z_xi, z_zeta
        real(rk) :: phi_xi, phi_eta, phi_zeta
        real(rk) :: q2, a2, q2_sonic

        ierr = 0

        x_xi   = node_deriv_xi(x, nj, nk, ni, i, j, k)
        y_xi   = node_deriv_xi(y, nj, nk, ni, i, j, k)
        z_xi   = node_deriv_xi(z, nj, nk, ni, i, j, k)
        y_eta  = node_deriv_eta(y, nj, nk, ni, i, j, k)
        z_zeta = node_deriv_zeta(z, nj, nk, ni, i, j, k)

        if (x_xi <= 0.0_rk .or. y_eta <= 0.0_rk .or. z_zeta <= 0.0_rk) then
            ierr = -1
            return
        end if

        phi_xi   = node_deriv_xi(phi, nj, nk, ni, i, j, k)
        phi_eta  = node_deriv_eta(phi, nj, nk, ni, i, j, k)
        phi_zeta = node_deriv_zeta(phi, nj, nk, ni, i, j, k)

        call physical_velocity(x_xi, y_xi, y_eta, z_xi, z_zeta, phi_xi, phi_eta, phi_zeta, u, v, w)

        q2_sonic = 2.0_rk / (gamma + 1.0_rk)
        q2 = u*u + v*v + w*w
        if (q2 >= q2_sonic) then
            ierr = -2
            return
        end if

        a2   = 1.0_rk - 0.5_rk * (gamma - 1.0_rk) * q2
        rho  = a2 ** (1.0_rk / (gamma - 1.0_rk))
        t    = a2
        p    = rho ** gamma
        mach = sqrt(q2 / a2)

    end subroutine node_output

    !> *************************************************************************
    !! * laplace_core_output (PUBLIC)
    !! *
    !! *   Node-centered physical output (README Sec8) at every node of the
    !! *   full domain (nj, nk, ni) -- no size reduction, unlike residual or
    !! *   metrics, since velocity/density/pressure/temperature/Mach are
    !! *   physically meaningful everywhere, boundaries included.
    !! *
    !! *   @param ni, nj, nk     Grid dimensions (streamwise, vertical, spanwise)
    !! *
    !! *   @param x, y, z        Node coordinates (nj, nk, ni)
    !! *
    !! *   @param phi            Converged potential field (nj, nk, ni)
    !! *
    !! *   @param gamma          Ratio of specific heats (> 1)
    !! *
    !! *   @param u, v, w        Nondimensional physical velocity (nj, nk, ni)
    !! *
    !! *   @param rho, p, t      Nondimensional density, pressure, temperature
    !! *                         (nj, nk, ni)
    !! *
    !! *   @param mach           Mach number (nj, nk, ni)
    !! *
    !! *   @param ierr           0 = OK; -1 = degenerate grid; -2 = sonic limit
    !! *                         reached; -3 = invalid gamma (<= 1)
    !! *
    !! ************************************************************************/
    pure subroutine laplace_core_output(ni, nj, nk, x, y, z, phi, gamma, u, v, w, rho, p, t, mach, ierr)
        integer(ik), intent(in)                         :: ni, nj, nk
        real(rk),    intent(in),  dimension(nj, nk, ni) :: x, y, z, phi
        real(rk),    intent(in)                         :: gamma
        real(rk),    intent(out), dimension(nj, nk, ni) :: u, v, w, rho, p, t, mach
        integer,     intent(out)                        :: ierr

        integer(ik) :: i, j, k
        integer     :: nerr

        ierr = 0

        if (gamma <= 1.0_rk) then
            ierr = -3
            return
        end if

        do i = 1, ni
            do k = 1, nk
                do j = 1, nj
                    call node_output(ni, nj, nk, x, y, z, phi, gamma, i, j, k, &
                                      u(j,k,i), v(j,k,i), w(j,k,i), &
                                      rho(j,k,i), p(j,k,i), t(j,k,i), mach(j,k,i), nerr)
                    if (nerr /= 0) then
                        ierr = nerr
                        return
                    end if
                end do
            end do
        end do

    end subroutine laplace_core_output

    pure subroutine thomas_solve(n, sub, diag, sup, rhs, x)
        !> (PRIVATE) Solves a tridiagonal system of n equations via the Thomas
        !! algorithm (forward elimination + back substitution). Sub(1) and sup(n)
        !! are not references -- no sub-diagonal entry on the first row, no
        !! super-diagonal entry on the last. Generic linear-algebra primitive.
        integer(ik), intent(in)                :: n
        real(rk),    intent(in),  dimension(n) :: sub, diag, sup, rhs
        real(rk),    intent(out), dimension(n) :: x

        real(rk), dimension(n) :: c_prime, d_prime
        real(rk)               :: m
        integer(ik)            :: i

        c_prime(1) = sup(1) / diag(1)
        d_prime(1) = rhs(1) / diag(1)

        do i = 2, n
            m          = diag(i) - sub(i) * c_prime(i-1)
            c_prime(i) = sup(i) / m
            d_prime(i) = (rhs(i) - sub(i) * d_prime(i-1)) / m
        end do

        x(n) = d_prime(n)
        do i = n - 1, 1, -1
            x(i) = d_prime(i) - c_prime(i) * x(i+1)
        end do

    end subroutine thomas_solve

    pure subroutine sweep(ni, nj, nk, phi, y, z, m_in, omega, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, max_r, rms_r)
        !> (PRIVATE) One SLOR sweep: for each (i, k) line, solves the tridiagonal
        !! defect-corection system a_S*dphi(j-1) - a_P*dphi(j) + a_N*dphi(j+1) = -R(i,j,k)
        !! -- derived by expanding the residual formula and substituting phi = phi_old + dphi;
        !! evertyhing except the pure eta-eta coupling cancels, leaving the residual itself
        !! as the RHS -- then over-relaxes: phi <- phi + omega*dphi. Mutates phi inplace.
        !!
        !! max_r/rms_r are a byproduct of the node_residual calls already made
        !! to build each line's RHS -- essentially free, but they describe phi
        !! as it was BEFORE this sweep's correction (and before any density
        !! update that follows it this iteration), not the just-updated phi.
        !! laplace_core_solve uses these for its per-iteration convergence
        !! check/progress print (a one-iteration lag is harmless there) and
        !! does one authoritative laplace_core_residual call after its loop
        !! exits, so the final max_r/rms_r it returns are exact.
        integer(ik), intent(in)                                 :: ni, nj, nk
        real(rk),    intent(inout), dimension(nj, nk, ni)       :: phi
        real(rk),    intent(in),    dimension(nj, nk, ni)       :: y, z
        real(rk),    intent(in)                                 :: m_in, omega
        real(rk),    intent(in),    dimension(nj,   nk,   ni-1) :: Aii, Aij, Aik
        real(rk),    intent(in),    dimension(nj-1, nk,   ni)   :: Aji, Ajj, Ajk
        real(rk),    intent(in),    dimension(nj,   nk-1, ni)   :: Aki, Akj, Akk
        real(rk),    intent(out)                                :: max_r, rms_r

        integer(ik) :: i, j, k
        real(rk)    :: a_s(nj), a_p(nj), a_n(nj), rhs(nj), dphi(nj)
        real(rk)    :: r_node, a_w, a_e, a_b, a_t
        real(rk)    :: sum_r2

        max_r  = 0.0_rk
        sum_r2 = 0.0_rk

        do i = 1, ni - 1
            do k = 1, nk
                do j = 1, nj

                    if (j > 1) then
                        a_s(j) = Ajj(j-1, k, i)
                    else
                        a_s(j) = 0.0_rk ! symmetry: no south neighbour, unused by thomas_solve anyway
                    end if

                    if (j < nj) then
                        a_n(j) = Ajj(j, k, i)
                    else
                        a_n(j) = 0.0_rk ! wall: no north neighbour, unused by thomas_solve anyway
                    end if

                    if (i > 1) then
                        a_w = Aii(j, k, i-1)
                    else
                        a_w = 0.0_rk ! inlet: prescribed constant flux, no phi-coefficient at all
                    end if

                    a_e = Aii(j, k, i) ! always valid: i+1 <= ni, even at the Dirichlet output plane

                    if (k > 1) then 
                        a_b = Akk(j, k-1, i)
                    else
                        a_b = 0.0_rk !symmetry
                    end if

                    if (k < nk) then
                        a_t = Akk(j, k, i)
                    else
                        a_t = 0.0_rk ! wall
                    end if

                    a_p(j) = a_w + a_e + a_s(j) + a_n(j) + a_b + a_t

                    call node_residual(ni, nj, nk, phi, y, z, m_in, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, i, j, k, r_node)

                    rhs(j) = -r_node
                    max_r  = max(max_r, abs(r_node))
                    sum_r2 = sum_r2 + r_node**2
                end do

                call thomas_solve(nj, a_s, -a_p, a_n, rhs, dphi)

                do j = 1, nj
                    phi(j, k, i) = phi(j, k, i) + omega * dphi(j)
                end do
            end do

        end do

        rms_r = sqrt(sum_r2 / real((ni - 1) * nk * nj, rk))

    end subroutine sweep

    pure subroutine laplace_core_solve(ni, nj, nk, x, y, z, phi, m_in, gamma, &
                                        omega, omega_rho, density_update_stride, &
                                        max_iter, tol, &
                                        n_iter_done, ierr, &
                                        progress_iter, max_r_history, rms_r_history, stop_flag)
        !> progress_iter/max_r_history/rms_r_history/stop_flag are all
        !! optional -- pass none of them (NULL from C) for an ordinary
        !! synchronous solve. All four exist for a caller running this on its
        !! own thread:
        !!
        !!   - progress_iter, if present, is updated via an atomic release-
        !!     store every iteration (src/laplace_atomic.c) -- a reader on
        !!     another thread should read it via an acquire-load (e.g.
        !!     aether's atomic_load_acq_u64); observing N there guarantees
        !!     max_r_history(1:N)/rms_r_history(1:N) are safe to read, even
        !!     if this routine has since raced ahead further.
        !!   - max_r_history/rms_r_history, if present, must be caller-
        !!     allocated arrays of length max_iter; slot i is sweep's (one-
        !!     iteration-behind, see sweep's own docs) residual norms after
        !!     iteration i.
        !!   - stop_flag, if present, is read every iteration via an atomic
        !!     acquire-load with the loop counter as a witness argument --
        !!     necessary so a PURE-licensed compiler can never treat repeated
        !!     calls as having "the same arguments" and cache/hoist the read
        !!     (verified in tests/fortran/test_atomic.f90). A caller sets it
        !!     nonzero (via a release-store) to request early exit; phi is
        !!     left at whatever sweep last produced, same as any other exit.
        integer(ik), intent(in)                                    :: ni, nj, nk
        real(rk),    intent(in),    dimension(nj, nk, ni)          :: x, y, z
        real(rk),    intent(inout), dimension(nj, nk, ni)          :: phi
        real(rk),    intent(in)                                    :: m_in, gamma
        real(rk),    intent(in)                                    :: omega, omega_rho
        integer(ik), intent(in)                                    :: density_update_stride
        integer(ik), intent(in)                                    :: max_iter
        real(rk),    intent(in)                                    :: tol
        integer(ik), intent(out)                                   :: n_iter_done
        integer,     intent(out)                                   :: ierr
        integer(ik), intent(inout), optional                       :: progress_iter
        real(rk),    intent(inout), optional, dimension(max_iter)  :: max_r_history, rms_r_history
        integer(ik), intent(in),    optional                       :: stop_flag

        integer(ik) :: iter
        integer     :: derr
        real(rk)    :: max_r, rms_r
        logical     :: stop_requested

        real(rk), allocatable, dimension(:,:,:) :: rho_xi,   rho_xi_new,   Aii, Aij, Aik, r
        real(rk), allocatable, dimension(:,:,:) :: rho_eta,  rho_eta_new,  Aji, Ajj, Ajk
        real(rk), allocatable, dimension(:,:,:) :: rho_zeta, rho_zeta_new, Aki, Akj, Akk

        ! raw grid-derivative terms (README Sec4): depend only on the fixed
        ! grid x,y,z, never on phi or rho, so computed exactly once here
        ! instead of on every density/metrics update inside the iteration loop
        real(rk), allocatable, dimension(:,:,:) :: xi_x_xi,   xi_y_xi,   xi_y_eta,   xi_z_xi,   xi_z_zeta
        real(rk), allocatable, dimension(:,:,:) :: eta_x_xi,  eta_y_xi,  eta_y_eta,  eta_z_xi,  eta_z_zeta
        real(rk), allocatable, dimension(:,:,:) :: zeta_x_xi, zeta_y_xi, zeta_y_eta, zeta_z_xi, zeta_z_zeta

        allocate(  rho_xi(nj,   nk, ni-1),   rho_xi_new(nj,   nk, ni-1), Aii(nj,   nk, ni-1), Aij(nj,   nk, ni-1), Aik(nj,   nk, ni-1) , r(nj, nk, ni-1) )
        allocate( rho_eta(nj-1, nk,   ni),  rho_eta_new(nj-1, nk,   ni), Aji(nj-1, nk,   ni), Ajj(nj-1, nk,   ni), Ajk(nj-1, nk,   ni) )
        allocate(rho_zeta(nj,   nk-1, ni), rho_zeta_new(nj,   nk-1, ni), Aki(nj,   nk-1, ni), Akj(nj,   nk-1, ni), Akk(nj,   nk-1, ni) )

        allocate(xi_x_xi(nj,nk,ni-1), xi_y_xi(nj,nk,ni-1), xi_y_eta(nj,nk,ni-1), xi_z_xi(nj,nk,ni-1), xi_z_zeta(nj,nk,ni-1))
        allocate(eta_x_xi(nj-1,nk,ni), eta_y_xi(nj-1,nk,ni), eta_y_eta(nj-1,nk,ni), eta_z_xi(nj-1,nk,ni), eta_z_zeta(nj-1,nk,ni))
        allocate(zeta_x_xi(nj,nk-1,ni), zeta_y_xi(nj,nk-1,ni), zeta_y_eta(nj,nk-1,ni), zeta_z_xi(nj,nk-1,ni), zeta_z_zeta(nj,nk-1,ni))

        ierr = 0
        n_iter_done = 0
        max_r = 0.0_rk
        rms_r = 0.0_rk

        ! grid-only raw metrics: computed once, reused for the whole solve
        call raw_metrics_xi  (ni, nj, nk, x, y, z, xi_x_xi,   xi_y_xi,   xi_y_eta,   xi_z_xi,   xi_z_zeta,   derr)
        if (derr == 0) call raw_metrics_eta (ni, nj, nk, x, y, z, eta_x_xi,  eta_y_xi,  eta_y_eta,  eta_z_xi,  eta_z_zeta,  derr)
        if (derr == 0) call raw_metrics_zeta(ni, nj, nk, x, y, z, zeta_x_xi, zeta_y_xi, zeta_y_eta, zeta_z_xi, zeta_z_zeta, derr)
        if (derr /= 0) then
            ierr = derr
            return
        end if

        ! initial density / metrics from the callers' starting guess for phi
        call density_xi_from_raw  (ni, nj, nk, phi, xi_x_xi,   xi_y_xi,   xi_y_eta,   xi_z_xi,   xi_z_zeta,   gamma, rho_xi,   derr)
        if (derr == 0) call density_eta_from_raw (ni, nj, nk, phi, eta_x_xi,  eta_y_xi,  eta_y_eta,  eta_z_xi,  eta_z_zeta,  gamma, rho_eta,  derr)
        if (derr == 0) call density_zeta_from_raw(ni, nj, nk, phi, zeta_x_xi, zeta_y_xi, zeta_y_eta, zeta_z_xi, zeta_z_zeta, gamma, rho_zeta, derr)
        if (derr /= 0) then
            ierr = derr
            return
        end if

        Aii =  rho_xi * xi_y_eta * xi_z_zeta / xi_x_xi
        Aij = -rho_xi * xi_y_xi  * xi_z_zeta / xi_x_xi
        Aik = -rho_xi * xi_z_xi  * xi_y_eta  / xi_x_xi

        Aji = -rho_eta * eta_y_xi * eta_z_zeta / eta_x_xi
        Ajj =  rho_eta * eta_x_xi * eta_z_zeta / eta_y_eta * (1.0_rk + (eta_y_xi / eta_x_xi)**2)
        Ajk =  rho_eta * eta_y_xi * eta_z_xi   / eta_x_xi

        Aki = -rho_zeta * zeta_z_xi * zeta_y_eta / zeta_x_xi
        Akj =  rho_zeta * zeta_y_xi * zeta_z_xi  / zeta_x_xi
        Akk =  rho_zeta * zeta_x_xi * zeta_y_eta / zeta_z_zeta * (1.0_rk + (zeta_z_xi / zeta_x_xi)**2)

        do iter = 1, max_iter
            call sweep(ni, nj, nk, phi, y, z, m_in, omega, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, max_r, rms_r)

            if (mod(iter, density_update_stride) == 0) then
                call density_xi_from_raw  (ni, nj, nk, phi, xi_x_xi,   xi_y_xi,   xi_y_eta,   xi_z_xi,   xi_z_zeta,   gamma, rho_xi_new,   derr)
                if (derr == 0) call density_eta_from_raw (ni, nj, nk, phi, eta_x_xi,  eta_y_xi,  eta_y_eta,  eta_z_xi,  eta_z_zeta,  gamma, rho_eta_new,  derr)
                if (derr == 0) call density_zeta_from_raw(ni, nj, nk, phi, zeta_x_xi, zeta_y_xi, zeta_y_eta, zeta_z_xi, zeta_z_zeta, gamma, rho_zeta_new, derr)
                if (derr /= 0) then
                    ierr = derr; n_iter_done = iter; return
                end if

                rho_xi   = rho_xi   + omega_rho * (rho_xi_new   - rho_xi)
                rho_eta  = rho_eta  + omega_rho * (rho_eta_new  - rho_eta)
                rho_zeta = rho_zeta + omega_rho * (rho_zeta_new - rho_zeta)

                Aii =  rho_xi * xi_y_eta * xi_z_zeta / xi_x_xi
                Aij = -rho_xi * xi_y_xi  * xi_z_zeta / xi_x_xi
                Aik = -rho_xi * xi_z_xi  * xi_y_eta  / xi_x_xi

                Aji = -rho_eta * eta_y_xi * eta_z_zeta / eta_x_xi
                Ajj =  rho_eta * eta_x_xi * eta_z_zeta / eta_y_eta * (1.0_rk + (eta_y_xi / eta_x_xi)**2)
                Ajk =  rho_eta * eta_y_xi * eta_z_xi   / eta_x_xi

                Aki = -rho_zeta * zeta_z_xi * zeta_y_eta / zeta_x_xi
                Akj =  rho_zeta * zeta_y_xi * zeta_z_xi  / zeta_x_xi
                Akk =  rho_zeta * zeta_x_xi * zeta_y_eta / zeta_z_zeta * (1.0_rk + (zeta_z_xi / zeta_x_xi)**2)

            end if

            n_iter_done = iter

            if (present(max_r_history)) max_r_history(iter) = max_r
            if (present(rms_r_history)) rms_r_history(iter) = rms_r
            if (present(progress_iter)) call laplace_atomic_store_rel_i64(progress_iter, iter)

            stop_requested = .false.
            if (present(stop_flag)) stop_requested = laplace_atomic_load_acq_i64(stop_flag, iter) /= 0_ik

            if (max_r < tol .or. stop_requested) exit

        end do

        ! sweep's max_r/rms_r (used above) are one iteration "behind" --
        ! cheap, since they're a byproduct of work sweep already does, but
        ! evaluated on phi/Aii..Akk as they were BEFORE the final iteration's
        ! correction (and density update), not what's actually being returned.
        ! One authoritative pass here, done once rather than every iteration,
        ! makes the final history entry (and progress_iter's last value) exact.
        call laplace_core_residual(ni, nj, nk, phi, y, z, m_in, Aii, Aij, Aik, Aji, Ajj, Ajk, Aki, Akj, Akk, r, derr)
        if (derr /= 0) then
            ierr = derr
            return
        end if

        max_r = maxval(abs(r))
        rms_r = sqrt(sum(r**2) / real(size(r), rk))

        if (n_iter_done > 0) then
            if (present(max_r_history)) max_r_history(n_iter_done) = max_r
            if (present(rms_r_history)) rms_r_history(n_iter_done) = rms_r
        end if
        if (present(progress_iter)) call laplace_atomic_store_rel_i64(progress_iter, n_iter_done)

    end subroutine laplace_core_solve


end module laplace_core

