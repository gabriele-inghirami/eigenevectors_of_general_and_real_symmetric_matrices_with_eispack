! This micro library computes the eigenvalues and the eigenvectors of
! general and real symmetric matrices.
!
! Each subroutine or function contains the description of its usage.
!
! This file contains the source code for the real symmetric case.
!
! These routines have been taken from https://netlib.org/eispack/ on Sept., 2022
! and minimally modified by Gabriele Inghirami (g.inghirami a-T gsi.de) to adapt
! their syntax from Fortran 77 to Fortran 90. However, there is no guarantee that
! this operation did not introduce bugs which were not present in the original
! version. Please, use the code exclusively at your own risk.
! Copyright: since the original version did not contain a Copyright note and it
! was developed by people working for public US institutions, it is reasonable
! to assume that this is PUBLIC DOMAIN code.

! .....................................................................

real(kind(1d0)) function pythag(a, b)
    
    implicit none
    real(kind(1d0)) :: a, b

! finds dsqrt(a**2+b**2) without overflow or destructive underflow

    double precision :: p, r, s, t, u
    
    p = dmax1(dabs(a),dabs(b))
    
    if (p .ne. 0.0d0) then
        r = (dmin1(dabs(a),dabs(b))/p)**2
        
        do while(.true.)
            t = 4.0d0 + r
            
            if (t .eq. 4.0d0) then
                exit
            end if
            
            s = r/t
            u = 1.0d0 + 2.0d0*s
            p = u*p
            r = (s/u)**2 * r
        end do
        
    end if
    
    pythag = p
    return

end function ! pythag

! .....................................................................

subroutine rs(nm, n, a, w, z, fv1, fv2, ierr)

    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer :: n, nm, ierr
    real(KIND=dp) :: a(nm,n), w(n), z(nm,n), fv1(n), fv2(n)

!   this subroutine calls the recommended sequence of
!   subroutines from the eigensystem subroutine package (eispack)
!   to find the eigenvalues and eigenvectors (if desired)
!   of a real symmetric matrix.
!
!   on input
!
!     nm  must be set to the row dimension of the two-dimensional
!     array parameters as declared in the calling program
!     dimension statement.
!
!     n  is the order of the matrix  a.
!
!     a  contains the real symmetric matrix.
!
!   on output
!
!     w  contains the eigenvalues in ascending order.
!
!     z  contains the eigenvectors if matz is not zero.
!
!     ierr  is an integer output variable set equal to an error
!     completion code described in the documentation for tqlrat
!     and tql2.  the normal completion code is zero.
!
!     fv1  and  fv2  are temporary storage arrays.
!
!   questions and comments should be directed to burton s. garbow,
!   mathematics and computer science div, argonne national laboratory
!
!   this version dated august 1983.
!
!   ------------------------------------------------------------------

    if (n .gt. nm) then
        ierr = 10 * n
        return 
    end if
    
!   .......... find both eigenvalues and eigenvectors ..........
    call  tred2(nm, n, a, w, fv1, z)
    call  tql2(nm, n, w, fv1, z, ierr)
      
    return

end subroutine ! rs

! .....................................................................
      
subroutine tql2(nm, n, d, e, z, ierr)

    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer :: i, j, k, l, m, n, ii, l1, l2, nm, mml, ierr
    real(KIND=dp) :: d(n), e(n), z(nm,n)
    real(KIND=dp) :: c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, tst1, tst2, pythag
    logical :: go_on

!   this subroutine is a translation of the algol procedure tql2,
!   num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!   wilkinson.
!   handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!   this subroutine finds the eigenvalues and eigenvectors
!   of a symmetric tridiagonal matrix by the ql method.
!   the eigenvectors of a full symmetric matrix can also
!   be found if  tred2  has been used to reduce this
!   full matrix to tridiagonal form.
!
!   on input
!
!     nm must be set to the row dimension of two-dimensional
!     array parameters as declared in the calling program
!     dimension statement.
!
!     n is the order of the matrix.
!
!     d contains the diagonal elements of the input matrix.
!
!     e contains the subdiagonal elements of the input matrix
!     in its last n-1 positions.  e(1) is arbitrary.
!
!     z contains the transformation matrix produced in the
!     reduction by  tred2, if performed.  if the eigenvectors
!     of the tridiagonal matrix are desired, z must contain
!     the identity matrix.
!
!   on output
!
!     d contains the eigenvalues in ascending order.  if an
!     error exit is made, the eigenvalues are correct but
!     unordered for indices 1,2,...,ierr-1.
!
!     e has been destroyed.
!
!     z contains orthonormal eigenvectors of the symmetric
!     tridiagonal (or full) matrix.  if an error exit is made,
!     z contains the eigenvectors associated with the stored
!     eigenvalues.
!
!     ierr is set to
!     zero       for normal return,
!     j          if the j-th eigenvalue has not been
!                determined after 30 iterations.
!
!   calls pythag for  dsqrt(a*a + b*b) .
!
!   questions and comments should be directed to burton s. garbow,
!   mathematics and computer science div, argonne national laboratory
!
!   this version dated august 1983.
!
!   ------------------------------------------------------------------

    ierr = 0
      
    if (n .eq. 1) then
        return
    end if

    do i = 2, n
          e(i-1) = e(i)
    end do

    f = 0.0d0
    tst1 = 0.0d0
    e(n) = 0.0d0

    do l = 1, n
        j = 0
        h = dabs(d(l)) + dabs(e(l))
         
        if (tst1 .lt. h) then
            tst1 = h
        end if
         
!   .......... look for small sub-diagonal element ..........
        do m = l, n
            tst2 = tst1 + dabs(e(m))
            
            if (tst2 .eq. tst1) then
                exit
            end if
            
!   .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
        end do ! 110

        if (m .ne. l) then
            go_on = .true.
             
            do while (go_on) 
             
                if (j .eq. 30) then
                 
!   .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
                    ierr = l
                    return
                end if
                 
                j = j + 1
                 
!   .......... form shift ..........
                l1 = l + 1
                l2 = l1 + 1
                g = d(l)
                p = (d(l1) - g) / (2.0d0 * e(l))
                r = pythag(p,1.0d0)
                d(l) = e(l) / (p + dsign(r,p))
                d(l1) = e(l) * (p + dsign(r,p))
                dl1 = d(l1)
                h = g - d(l)
                 
                if (l2 .le. n) then
                
                    do i = l2, n
                        d(i) = d(i) - h
                    end do
                     
                end if
    
                f = f + h
                 
!   .......... ql transformation ..........
                p = d(m)
                c = 1.0d0
                c2 = c
                el1 = e(l1)
                s = 0.0d0
                mml = m - l
                 
!   .......... for i=m-1 step -1 until l do -- ..........
                do ii = 1, mml
                    c3 = c2
                    c2 = c
                    s2 = s
                    i = m - ii
                    g = c * e(i)
                    h = c * p
                    r = pythag(p,e(i))
                    e(i+1) = s * r
                    s = e(i) / r
                    c = p / r
                    p = c * d(i) - s * g
                    d(i+1) = h + s * (c * g + s * d(i))
                    
!   .......... form vector ..........
                    do k = 1, n
                        h = z(k,i+1)
                        z(k,i+1) = s * z(k,i) + c * h
                        z(k,i) = c * z(k,i) - s * h
                    end do
                    
                end do

                p = -s * s2 * c3 * el1 * e(l) / dl1
                e(l) = s * p
                d(l) = c * p
                tst2 = tst1 + dabs(e(l))
             
                if (tst2 .le. tst1) then
                    go_on = .false.
                end if
             
            end do ! end go on
        
        end if
     
        d(l) = d(l) + f
    end do

!   .......... order eigenvalues and eigenvectors ..........
    do ii = 2, n
        i = ii - 1
        k = i
        p = d(i)

        do j = ii, n
        
            if (d(j) .ge. p) then
                cycle
            end if
            
            k = j
            p = d(j)
        end do

        if (k .eq. i) then
            cycle
        end if 
        
        d(k) = d(i)
        d(i) = p

        do j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
        end do
        
    end do

    return
    
end subroutine ! tql2
      
! .....................................................................

subroutine tred2(nm, n, a, d, e, z)

    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer :: i, j, k, l, n, ii, nm, jp1
    real(KIND=dp) :: a(nm,n), d(n), e(n), z(nm,n)
    real(KIND=dp) :: f, g, h, hh, scaling

!   this subroutine is a translation of the algol procedure tred2,
!   num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!   handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!   this subroutine reduces a real symmetric matrix to a
!   symmetric tridiagonal matrix using and accumulating
!   orthogonal similarity transformations.
!
!   on input
!
!     nm must be set to the row dimension of two-dimensional
!     array parameters as declared in the calling program
!     dimension statement.
!
!     n is the order of the matrix.
!
!     a contains the real symmetric input matrix.  only the
!     lower triangle of the matrix need be supplied.
!
!   on output
!
!     d contains the diagonal elements of the tridiagonal matrix.
!
!     e contains the subdiagonal elements of the tridiagonal
!     matrix in its last n-1 positions.  e(1) is set to zero.
!
!     z contains the orthogonal transformation matrix
!     produced in the reduction.
!
!     a and z may coincide.  if distinct, a is unaltered.
!
!   questions and comments should be directed to burton s. garbow,
!   mathematics and computer science div, argonne national laboratory
!
!   this version dated august 1983.
!
!   ------------------------------------------------------------------

    do i = 1, n
    
        do j = i, n
            z(j,i) = a(j,i)
        end do
        
        d(i) = a(n,i)
    end do

    if (n .ne. 1) then
    
!   .......... for i=n step -1 until 2 do -- ..........
        do ii = 2, n
            i = n + 2 - ii
            l = i - 1
            h = 0.0d0
            scaling = 0.0d0
            
            if (l .ge. 2) then
            
!     .......... scale row (algol tol then not needed) ..........
                do k = 1, l
                    scaling = scaling + dabs(d(k))
                end do
                
            end if
            
            if ((( l .ge. 2) .and. (scaling .eq. 0.0d0)) .or. (l .lt. 2)) then
                e(i) = d(l)
                do j = 1, l
                    d(j) = z(l,j)
                    z(i,j) = 0.0d0
                    z(j,i) = 0.0d0
                end do

            else

                do k = 1, l
                    d(k) = d(k) / scaling
                    h = h + d(k) * d(k)
                end do
  
                f = d(l)
                g = -dsign(dsqrt(h),f)
                e(i) = scaling * g
                h = h - f * g
                d(l) = f - g
            
!   .......... form a*u ..........
                do j = 1, l 
                    e(j) = 0.0d0
                end do

                do j = 1, l
                    f = d(j)
                    z(j,i) = f
                    g = e(j) + z(j,j) * f
                    jp1 = j + 1
                
                    if (l .ge. jp1) then
                
                        do k = jp1, l
                            g = g + z(k,j) * d(k)
                            e(k) = e(k) + z(k,j) * f
                        end do
                    
                    end if 
                
                    e(j) = g
                end do
            
!   .......... form p ..........
                f = 0.0d0

                do j = 1, l
                    e(j) = e(j) / h
                    f = f + e(j) * d(j)
                end do

                hh = f / (h + h)
            
!   .......... form q ..........
                do j = 1, l
                    e(j) = e(j) - hh * d(j)
                end do
            
!   .......... form reduced a ..........
                do j = 1, l
                    f = d(j)
                    g = e(j)
                
                    do k = j, l
                        z(k,j) = z(k,j) - f * e(k) - g * d(k)
                    end do
                
                    d(j) = z(l,j)
                    z(i,j) = 0.0d0
                end do
           
            end if
        
            d(i) = h
        end do

!   .......... accumulation of transformation matrices ..........
        do i = 2, n
            l = i - 1
            z(n,l) = z(l,l)
            z(l,l) = 1.0d0
            h = d(i)
        
            if (h .ne. 0.0d0) then
        
                do k = 1, l
                    d(k) = z(k,i) / h
                end do

                do j = 1, l
                    g = 0.0d0
                
                    do k = 1, l
                        g = g + z(k,i) * z(k,j)
                    end do
                
                    do k = 1, l
                        z(k,j) = z(k,j) - g * d(k)
                    end do
                
                end do
            
            end if
        
            do k = 1, l
                z(k,i) = 0.0d0
            end do

        end do
    end if 
    do i = 1, n
        d(i) = z(n,i)
        z(n,i) = 0.0d0
    end do

    z(n,n) = 1.0d0
    e(1) = 0.0d0
    return
end subroutine ! tred2
