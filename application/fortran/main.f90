!> Native reference driver. Calls only laplace_f_api, never laplace_core
!! directly, so this also proves the f_api surface is sufficient to drive
!! the whole program.
program laplace_fortran_app
    use, intrinsic :: iso_fortran_env, only : output_unit
    use :: laplace_types,  only : rk, ik
    use :: laplace_f_api,  only : laplace_kernel

    implicit none

    integer(ik), parameter :: n = 5_ik
    real(rk) :: x(n), y(n)
    integer(ik) :: i

    do i = 1, n
        x(i) = real(i, rk)
    end do

    call laplace_kernel(n, x, y)

    write(output_unit, '(A)') 'laplace fortran_app: placeholder kernel smoke test'
    do i = 1, n
        write(output_unit, '(I3, 2F10.4)') i, x(i), y(i)
    end do

end program laplace_fortran_app
