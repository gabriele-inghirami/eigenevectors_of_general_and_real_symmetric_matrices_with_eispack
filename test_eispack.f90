! Code to test the mini library, extracted from EISPACK, composed by the source files:
! eispack_symm_matrix.f90 and eispack_general_matrix.f90.
!
! For an overview of what this example/test code does, please, read the file README.md.
! 
! Author: Gabriele Inghirami ( gabriele.inghirami a.  gmail.com)
! License: PUBLIC DOMAIN

subroutine display_eigenvalues_and_eigenvectors(n, eva, eve)
    
    implicit none
    integer :: i, j
    integer, intent(in) :: n
    real(kind(1.0d0)), dimension(n), intent(in) :: eva
    real(kind(1.0d0)), dimension(n,n), intent(in) :: eve
    character(len=*), parameter  :: fmt1 = "(a11,i2,a1,f15.10)"
    character(len=*), parameter  :: fmt2 = "(a20,i2,a1,4f15.10)"

    write(*,*)
    do i = 1, n
        write(*,fmt1,advance="no") " Eigenvalue", i, ":", eva(i)
        write(*,fmt2,advance="no") "Eigenvector",i,":", (eve(j,i), j = 1, n)
        write(*,*)
    end do
    write(*,*)

end subroutine

subroutine display_square_matrix(n, mat)
    
    implicit none
    integer :: i, j
    integer, intent(in) :: n
    real(kind(1.0d0)), dimension(1:n,1:n), intent(in) :: mat 
    character(len=*), parameter  :: fmt1 = "(4f15.10)"

    write(*,*)
    do i = 1, n
        write(*,fmt1,advance="no") (mat(i,j), j = 1, n)
        write(*,*)
    end do
    write(*,*)

end subroutine

subroutine square_matrix_vector_dot_product(mat, vec, n, res)
    implicit none
    integer, intent(in) :: n
    real(kind(1.0d0)), dimension(1:n,1:n), intent(in) :: mat
    real(kind(1.0d0)), dimension(1:n), intent(in) :: vec
    real(kind(1.0d0)), dimension(1:n), intent(out) :: res
    integer :: i, j
    res(:) = 0.d0
    do i = 1, n
        do j = 1, n
            res(i) = res(i) + mat(i,j)*vec(j)
        end do
    end do
end subroutine

subroutine print_symmetric_mat_case_check_message
    write(*,*) "Check of the solution."
    write(*,*) "Maximum absolute value among the components of ( T . x - lambda x ) for eigenvalue lambda and &
                eigenvector x:"//NEW_LINE('a')
end subroutine

subroutine print_generalised_case_check_message
    write(*,*) "Check of the solution."
    write(*,*) "Maximum absolute value among the components of ( T . x - lambda B . x ) for eigenvalue lambda and &
                eigenvector x:"//NEW_LINE('a')
end subroutine

program test
    
    implicit none
    integer, parameter :: n = 4
    integer :: i, j, ierr
    real(kind(1.0d0)) :: lambda
    real(kind(1.0d0)), dimension(n) :: w, fv1, fv2, alfr, alfi, beta, res_left, res_right, res_final, vec
    real(kind(1.0d0)), dimension(n,n) :: a, b, z, a_copy, b_copy

    write(*,*) NEW_LINE('a')//"Test 1"//NEW_LINE('a')
    write(*,*) "Testing real symmetrix matrix eigenvalues and eigenvectors solver with:"
    ! values are inserted in the matrix in the order a11, a21, a31, ..., an1, a12, a22, a32,...
    a = reshape((/3d0,2d0,-1d0,6d0,2d0,0.5d0,0.4d0,1.7d0,-1d0,0.4d0,1.9d0,4.1d0,6d0,1.7d0,4.1d0,1.5d0/),shape(a))
    a_copy = a
    call display_square_matrix(n, a)
    call rs(n, n, a, w, z, fv1, fv2, ierr)
    write(*,*) "Solution:"
    call display_eigenvalues_and_eigenvectors(n, w, z)
    do i = 1, n
        do j = 1, n
            vec(j) = z(j,i)
        end do
        lambda = w(i)
        call square_matrix_vector_dot_product(a_copy, vec, n, res_left)
        res_final = res_left - lambda * vec
        write(*,*) i," : ", maxval(abs(res_final))
    end do

    write(*,*) NEW_LINE('a')//"Test 2"//NEW_LINE('a')
    write(*,*) "Testing real symmetrix matrix eigenvalues and eigenvectors solver with:"
    ! values are inserted in the matrix in the order a11, a21, a31, ..., an1, a12, a22, a32,...
    a = reshape((/21d0,-3.8d0,2.6d0,5d0,-3.8d0,4.4d0,-8.5d0,-48d0,2.6d0,-8.5d0,5.9d0,3.2d0,5d0,-48d0,3.2d0,9.9d0/),shape(a))
    a_copy = a
    call display_square_matrix(n, a)
    call rs(n, n, a, w, z, fv1, fv2, ierr)
    write(*,*) "Solution:"
    call display_eigenvalues_and_eigenvectors(n, w, z)
    call print_symmetric_mat_case_check_message
    do i = 1, n
        do j = 1, n
            vec(j) = z(j,i)
        end do
        lambda = w(i)
        call square_matrix_vector_dot_product(a_copy, vec, n, res_left)
        res_final = res_left - lambda * vec
        write(*,*) i," : ", maxval(abs(res_final))
    end do

    write(*,*) NEW_LINE('a')//"Test 3"//NEW_LINE('a')
    write(*,*) "Testing generalized eigenvalues and eigenvectors solver (fot the eq. A . x = lambda B . x) with A and B:"
    ! values are inserted in the matrix in the order a11, a21, a31, ..., an1, a12, a22, a32,...
    a = reshape((/3d0,2d0,-1d0,6d0,2d0,0.5d0,0.4d0,1.7d0,-1d0,0.4d0,1.9d0,4.1d0,6d0,1.7d0,4.1d0,1.5d0/),shape(a))
    a = reshape((/1.1d0,0.23d0,0.68d0,0.295d0,0.2d0,-0.14d0,0.39d0,0.56d0,0.7d0,0.41d0,1.03d0,0.91d0,&
        0.34d0,0.51d0,0.87d0,1.6d0/),shape(a))
    b = reshape((/1.d0,0d0,0d0,0d0,0d0,1.d0,0d0,0d0,0d0,0d0,1.d0,0d0,0d0,0d0,0d0,1.d0/),shape(b))
    ! we make copies that will be used at the end to check the results because the original matrices will be modified
    a_copy = a
    b_copy = b
    call display_square_matrix(n, b)
    call display_square_matrix(n, a)
    call rgg(n, n, a, b, alfr, alfi, beta, 1, z, ierr)
    if (ierr .ne. 0) then
        write(*,*) "Error in Eispack subroutine rgg\n"
        call exit(1)
    end if
    write(*,*) "Solution:"
    ! warning: here we display the solution assuming that the values are all real
    call display_eigenvalues_and_eigenvectors(n, alfr/beta, z)
    call print_generalised_case_check_message
    do i = 1, n
        do j = 1, n
            vec(j) = z(j,i)
        end do
        lambda = alfr(i)/beta(i)
        call square_matrix_vector_dot_product(a_copy, vec, n, res_left)
        call square_matrix_vector_dot_product(b_copy, vec, n, res_right)
        res_final = res_left - lambda * res_right
        write(*,*) i," : ", maxval(abs(res_final))
    end do
            
    write(*,*) NEW_LINE('a')//"Test 4"//NEW_LINE('a')
    write(*,*) "Testing generalized eigenvalues and eigenvectors solver for the eq. ( A . x = lambda B . x ) with A and B:"
    ! values are inserted in the matrix in the order a11, a21, a31, ..., an1, a12, a22, a32,...
    a = reshape((/3d0,2d0,-1d0,6d0,2d0,0.5d0,0.4d0,1.7d0,-1d0,0.4d0,1.9d0,4.1d0,6d0,1.7d0,4.1d0,1.5d0/),shape(a))
    a = reshape((/3.2d0,0.1d0,0.4d0,0.5d0,0.13d0,-0.25d0,0.22d0,0.91d0,0.38d0,0.21d0,0.2d0,0.08d0,&
        0.47d0,0.9d0,0.1d0,1.4d0/),shape(a))
    b = reshape((/1.d0,0d0,0d0,0d0,0d0,1.d0,0d0,0d0,0d0,0d0,1.d0,0d0,0d0,0d0,0d0,1.d0/),shape(b))
    ! we make copies that will be used at the end to check the results because the original matrices will be modified
    a_copy = a
    b_copy = b
    call display_square_matrix(n, b)
    call display_square_matrix(n, a)
    call rgg(n, n, a, b, alfr, alfi, beta, 1, z, ierr)
    if (ierr .ne. 0) then
        write(*,*) "Error in Eispack subroutine rgg\n"
        call exit(1)
    end if
    write(*,*) "Solution:"
    ! warning: here we display the solution assuming that the values are all real
    call display_eigenvalues_and_eigenvectors(n, alfr/beta, z)
    call print_generalised_case_check_message
    do i = 1, n
        do j = 1, n
            vec(j) = z(j,i)
        end do 
        lambda = alfr(i)/beta(i)
        call square_matrix_vector_dot_product(a_copy, vec, n, res_left)
        call square_matrix_vector_dot_product(b_copy, vec, n, res_right)
        res_final = res_left - lambda * res_right
        write(*,*) i," : ", maxval(abs(res_final))
    end do

    write(*,*) NEW_LINE('a')//"Test 5"//NEW_LINE('a')
    write(*,*) "Testing generalized eigenvalues and eigenvectors solver (fot the eq. A . x = lambda B . x) with A and B:"
    ! values are inserted in the matrix in the order a11, a21, a31, ..., an1, a12, a22, a32,...
    a = reshape((/3d0,2d0,-1d0,6d0,2d0,0.5d0,0.4d0,1.7d0,-1d0,0.4d0,1.9d0,4.1d0,6d0,1.7d0,4.1d0,1.5d0/),shape(a))
    a = reshape((/-5d0,1d0,-1.5d0,2.7d0,1d0,2d0,4.7d0,5.1d0,-1d0,4d0,0.4d0,3.3d0,2.6d0,1.1d0,0.5d0,-1d0/),shape(a))
    b = reshape((/1.2d0,0.11d0,0.05d0,0.03d0,0.1d0,0.95d0,0.09d0,0.08d0,0.2d0,0.25d0,1.15d0,0.075d0,&
                  0.15d0,0.12d0,0.14d0,1.05d0/),shape(b))
    ! we make copies that will be used at the end to check the results because the original matrices will be modified
    a_copy = a
    b_copy = b
    call display_square_matrix(n, b)
    call display_square_matrix(n, a)
    call rgg(n, n, a, b, alfr, alfi, beta, 1, z, ierr)
    if (ierr .ne. 0) then
        write(*,*) "Error in Eispack subroutine rgg\n"
        call exit(1)
    end if
    write(*,*) "Solution:"
    ! warning: here we display the solution assuming that the values are all real
    call display_eigenvalues_and_eigenvectors(n, alfr/beta, z)
    call print_generalised_case_check_message
    do i = 1, n
        do j = 1, n
            vec(j) = z(j,i)
        end do 
        lambda = alfr(i)/beta(i)
        call square_matrix_vector_dot_product(a_copy, vec, n, res_left)
        call square_matrix_vector_dot_product(b_copy, vec, n, res_right)
        res_final = res_left - lambda * res_right
        write(*,*) i," : ", maxval(abs(res_final))
    end do
    write(*,*) ""

end program

