module laplace_types
    use, intrinsic :: iso_c_binding, only : c_double, c_int64_t

    implicit none

    private

    public :: rk
    public :: ik

    integer, parameter :: rk = c_double
    integer, parameter :: ik = c_int64_t


end module laplace_types
