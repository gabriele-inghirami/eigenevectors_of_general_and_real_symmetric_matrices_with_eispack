! Code to test the mini library, extracted from EISPACK, eispack_symm_matrix.f90.
! Author: Gabriele Inghirami ( gabriele.inghirami a.  gmail.com)
! License: PUBLIC DOMAIN

subroutine display(a, b, eva, eve)
    
    implicit none
    integer :: i, j
    integer, intent(in) :: a, b
    real(kind(1.0d0)), dimension(a), intent(in) :: eva
    real(kind(1.0d0)), dimension(a,b), intent(in) :: eve
    character(len=*), parameter  :: fmt1 = "(a10,i2,a1,f14.9)"
    character(len=*), parameter  :: fmt2 = "(a17,4f14.9)"

    write(*,*)
    do i = 1, a
        write(*,fmt1,advance="no") "Eigenvalue", i, ":", eva(i)
        write(*,fmt2,advance="no") "Eigenvector:", (eve(j,i), j = 1, b)
        write(*,*)
    end do
    write(*,*)

end subroutine

program test
    
    implicit none
    integer, parameter :: n = 4, nm = n
    integer :: ierr
    real(kind(1.0d0)), dimension(n) :: w, fv1, fv2
    real(kind(1.0d0)), dimension(nm,n) :: a(nm,n), z(nm,n)

    a = reshape((/3d0,2d0,-1d0,6d0,2d0,0.5d0,0.4d0,1.7d0,-1d0,0.4d0,1.9d0,4.1d0,6d0,1.7d0,4.1d0,1.5d0/),shape(a))

    call rs(nm, n, a, w, z, fv1, fv2, ierr)
    call display(nm, n, w, z)

    a = reshape((/21d0,-3.8d0,2.6d0,5d0,-3.8d0,4.4d0,-8.5d0,-48d0,2.6d0,-8.5d0,5.9d0,3.2d0,5d0,-48d0,3.2d0,9.9d0/),shape(a))

    call rs(nm, n, a, w, z, fv1, fv2, ierr)
    call display(nm, n, w, z)

end program

