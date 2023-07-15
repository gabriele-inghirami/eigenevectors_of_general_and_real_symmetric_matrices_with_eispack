! This tiny library computes the eigenvalues and the eigenvectors of
! real symmetric matrices and for the generalized eigenvalue problem
! A . x = lambda B . x, where A and B are general matrices.
!
! Each subroutine or function contains the description of its usage.
!
! This file contains the source code for the .
!
! These routines have been taken from https://netlib.org/eispack/ on Sept., 2022
! and minimally modified by Gabriele Inghirami (g.inghirami a-T gsi.de) to adapt
! their syntax from Fortran 77 to Fortran 90. However, there is no guarantee that
! this operation did not introduce bugs which were not present in the original
! version. Please, use the code exclusively at your own risk.
! Copyright: since the original version did not contain a Copyright note and it
! was developed by people working for public US institutions, it is reasonable
! to assume that this is PUBLIC DOMAIN code.
!
! .....................................................................

subroutine qzhes(nm,n,a,b,matz,z)
    implicit none

    integer i,j,k,l,n,lb,l1,nm,nk1,nm1,nm2
    real(kind(1.0d0)) a(nm,n),b(nm,n),z(nm,n)
    real(kind(1.0d0)) r,s,t,u1,u2,v1,v2,rho
    logical matz

!   this subroutine is the first step of the qz algorithm
!   for solving generalized matrix eigenvalue problems,
!   siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!   this subroutine accepts a pair of real general matrices and
!   reduces one of them to upper hessenberg form and the other
!   to upper triangular form using orthogonal transformations.
!   it is usually followed by  qzit,  qzval  and, possibly,  qzvec.
!
!   on input
!
!     nm must be set to the row dimension of two-dimensional
!     array parameters as declared in the calling program
!     dimension statement.
!
!     n is the order of the matrices.
!
!     a contains a real general matrix.
!
!     b contains a real general matrix.
!
!     matz should be set to .true. if the right hand transformations
!     are to be accumulated for later use in computing
!     eigenvectors, and to .false. otherwise.
!
!     on output
!
!     a has been reduced to upper hessenberg form.  the elements
!     below the first subdiagonal have been set to zero.
!
!     b has been reduced to upper triangular form.  the elements
!     below the main diagonal have been set to zero.
!
!     z contains the product of the right hand transformations if
!     matz has been set to .true.  otherwise, z is not referenced.
!
!   questions and comments should be directed to burton s. garbow,
!   mathematics and computer science div, argonne national laboratory
!
!   this version dated august 1983.
!
!     ------------------------------------------------------------------

!   .......... initialize z ..........
    if (matz) then
        do j = 1, n
            do i = 1, n
                z(i,j) = 0.0d0
            end do
            z(j,j) = 1.0d0
        end do
    end if
!   .......... reduce b to upper triangular form ..........
    if (n .le. 1) return
    nm1 = n - 1
    do l = 1, nm1
        l1 = l + 1
        s = 0.0d0
        do i = l1, n
            s = s + dabs(b(i,l))
        end do
        if (s .eq. 0.0d0) cycle
        s = s + dabs(b(l,l))
        r = 0.0d0
        do i = l, n
            b(i,l) = b(i,l) / s
            r = r + b(i,l)**2
        end do
        r = dsign(dsqrt(r),b(l,l))
        b(l,l) = b(l,l) + r
        rho = r * b(l,l)
        do j = l1, n
            t = 0.0d0
            do i = l, n
                t = t + b(i,l) * b(i,j)
            end do
            t = -t / rho
            do i = l, n
                b(i,j) = b(i,j) + t * b(i,l)
            end do
        end do
        do j = 1, n
            t = 0.0d0
            do i = l, n
                t = t + b(i,l) * a(i,j)
            end do
            t = -t / rho
            do i = l, n
                a(i,j) = a(i,j) + t * b(i,l)
            end do
        end do
        b(l,l) = -s * r
        do i = l1, n
            b(i,l) = 0.0d0
        end do
    end do
!   .......... reduce a to upper hessenberg form, while
!              keeping b triangular ..........
    if (n .eq. 2) return
    nm2 = n - 2
    do k = 1, nm2
        nk1 = nm1 - k
!   .......... for l=n-1 step -1 until k+1 do -- ..........
        do lb = 1, nk1
            l = n - lb
            l1 = l + 1
!   .......... zero a(l+1,k) ..........
            s = dabs(a(l,k)) + dabs(a(l1,k))
            if (s .eq. 0.0d0) cycle
            u1 = a(l,k) / s
            u2 = a(l1,k) / s
            r = dsign(dsqrt(u1*u1+u2*u2),u1)
            v1 =  -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1
            do j = k, n
                t = a(l,j) + u2 * a(l1,j)
                a(l,j) = a(l,j) + t * v1
                a(l1,j) = a(l1,j) + t * v2
            end do
            a(l1,k) = 0.0d0
            do j = l, n
                t = b(l,j) + u2 * b(l1,j)
                b(l,j) = b(l,j) + t * v1
                b(l1,j) = b(l1,j) + t * v2
            end do
!   .......... zero b(l+1,l) ..........
            s = dabs(b(l1,l1)) + dabs(b(l1,l))
            if (s .eq. 0.0d0) cycle
            u1 = b(l1,l1) / s
            u2 = b(l1,l) / s
            r = dsign(dsqrt(u1*u1+u2*u2),u1)
            v1 =  -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1
            do i = 1, l1
                t = b(i,l1) + u2 * b(i,l)
                b(i,l1) = b(i,l1) + t * v1
                b(i,l) = b(i,l) + t * v2
            end do
            b(l1,l) = 0.0d0
            do i = 1, n
                t = a(i,l1) + u2 * a(i,l)
                a(i,l1) = a(i,l1) + t * v1
                a(i,l) = a(i,l) + t * v2
            end do
            if (.not. matz) cycle
            do i = 1, n
                t = z(i,l1) + u2 * z(i,l)
                z(i,l1) = z(i,l1) + t * v1
                z(i,l) = z(i,l) + t * v2
            end do
        end do
    end do
end subroutine


subroutine qzit(nm,n,a,b,eps1,matz,z,ierr)
    implicit none

    integer i,j,k,l,n,en,k1,k2,ld,ll,l1,na,nm,ish,itn,its,km1,lm1
    integer enm2,ierr,lor1,enorn
    real(kind(1.0d0)) a(nm,n),b(nm,n),z(nm,n)
    real(kind(1.0d0)) r,s,t,a1,a2,a3,ep,sh,u1,u2,u3,v1,v2,v3,ani,a11
    real(kind(1.0d0)) a12,a21,a22,a33,a34,a43,a44,bni,b11,b12,b22,b33,b34
    real(kind(1.0d0)) b44,epsa,epsb,eps1,anorm,bnorm
    logical matz,notlas

!   this subroutine is the second step of the qz algorithm
!   for solving generalized matrix eigenvalue problems,
!   siam j. numer. anal. 10, 241-256(1973) by moler and stewart,
!   as modified in technical note nasa tn d-7305(1973) by ward.
!
!   this subroutine accepts a pair of real matrices, one of them
!   in upper hessenberg form and the other in upper triangular form.
!   it reduces the hessenberg matrix to quasi-triangular form using
!   orthogonal transformations while maintaining the triangular form
!   of the other matrix.  it is usually preceded by  qzhes  and
!   followed by  qzval  and, possibly,  qzvec.
!
!   on input
!
!     nm must be set to the row dimension of two-dimensional
!     array parameters as declared in the calling program
!     dimension statement.
!
!     n is the order of the matrices.
!
!     a contains a real upper hessenberg matrix.
!
!     b contains a real upper triangular matrix.
!
!     eps1 is a tolerance used to determine negligible elements.
!     eps1 = 0.0 (or negative) may be input, in which case an
!     element will be neglected only if it is less than roundoff
!     error times the norm of its matrix.  if the input eps1 is
!     positive, then an element will be considered negligible
!     if it is less than eps1 times the norm of its matrix.  a
!     positive value of eps1 may result in faster execution,
!     but less accurate results.
!
!     matz should be set to .true. if the right hand transformations
!     are to be accumulated for later use in computing
!     eigenvectors, and to .false. otherwise.
!
!     z contains, if matz has been set to .true., the
!     transformation matrix produced in the reduction
!     by  qzhes, if performed, or else the identity matrix.
!     if matz has been set to .false., z is not referenced.
!
!   on output
!
!     a has been reduced to quasi-triangular form.  the elements
!     below the first subdiagonal are still zero and no two
!     consecutive subdiagonal elements are nonzero.
!
!     b is still in upper triangular form, although its elements
!     have been altered.  the location b(n,1) is used to store
!     eps1 times the norm of b for later use by  qzval  and  qzvec.
!
!     z contains the product of the right hand transformations
!     (for both steps) if matz has been set to .true..
!
!     ierr is set to
!       zero       for normal return,
!       j          if the limit of 30*n iterations is exhausted
!                while the j-th eigenvalue is being sought.
!
!   questions and comments should be directed to burton s. garbow,
!   mathematics and computer science div, argonne national laboratory
!
!   this version dated august 1983.
!
!   ------------------------------------------------------------------

    ierr = 0
!   .......... compute epsa,epsb ..........
    anorm = 0.0d0
    bnorm = 0.0d0

    do i = 1, n
        ani = 0.0d0
        if (i .ne. 1) ani = dabs(a(i,i-1))
        bni = 0.0d0
        do j = i, n
            ani = ani + dabs(a(i,j))
            bni = bni + dabs(b(i,j))
        end do
        if (ani .gt. anorm) anorm = ani
        if (bni .gt. bnorm) bnorm = bni
    end do
    if (anorm .eq. 0.0d0) anorm = 1.0d0
    if (bnorm .eq. 0.0d0) bnorm = 1.0d0
    ep = eps1
    if (ep .gt. 0.0d0) goto 50
!   .......... use roundoff level if eps1 is zero ..........
    ep = epsilon(1.0d0)
 50 epsa = ep * anorm
    epsb = ep * bnorm
!   .......... reduce a to quasi-triangular form, while
!              keeping b triangular ..........
    lor1 = 1
    enorn = n
    en = n
    itn = 30*n
!   .......... begin qz step ..........
 60 if (en .le. 2) goto 999
    if (.not. matz) enorn = en
    its = 0
    na = en - 1
    enm2 = na - 1
 70 ish = 2
!   .......... check for convergence or reducibility.
!              for l=en step -1 until 1 do -- ..........
    do ll = 1, en
        lm1 = en - ll
        l = lm1 + 1
        if (l .eq. 1) goto 95
        if (dabs(a(l,lm1)) .le. epsa) goto 90
    end do
 90 a(l,lm1) = 0.0d0
    if (l .lt. na) goto 95
!   .......... 1-by-1 or 2-by-2 block isolated ..........
    en = lm1
    goto 60
!   .......... check for small top of b ..........
 95 ld = l
100 l1 = l + 1
    b11 = b(l,l)
    if (dabs(b11) .gt. epsb) goto 120
    b(l,l) = 0.0d0
    s = dabs(a(l,l)) + dabs(a(l1,l))
    u1 = a(l,l) / s
    u2 = a(l1,l) / s
    r = dsign(dsqrt(u1*u1+u2*u2),u1)
    v1 = -(u1 + r) / r
    v2 = -u2 / r
    u2 = v2 / v1
    do j = l, enorn
        t = a(l,j) + u2 * a(l1,j)
        a(l,j) = a(l,j) + t * v1
        a(l1,j) = a(l1,j) + t * v2
        t = b(l,j) + u2 * b(l1,j)
        b(l,j) = b(l,j) + t * v1
        b(l1,j) = b(l1,j) + t * v2
    end do
    if (l .ne. 1) a(l,lm1) = -a(l,lm1)
    lm1 = l
    l = l1
    goto 90
120 a11 = a(l,l) / b11
    a21 = a(l1,l) / b11
    if (ish .eq. 1) goto 140
!   .......... iteration strategy ..........
    if (itn .eq. 0) goto 998
    if (its .eq. 10) goto 155
!   .......... determine type of shift ..........
    b22 = b(l1,l1)
    if (dabs(b22) .lt. epsb) b22 = epsb
    b33 = b(na,na)
    if (dabs(b33) .lt. epsb) b33 = epsb
    b44 = b(en,en)
    if (dabs(b44) .lt. epsb) b44 = epsb
    a33 = a(na,na) / b33
    a34 = a(na,en) / b44
    a43 = a(en,na) / b33
    a44 = a(en,en) / b44
    b34 = b(na,en) / b44
    t = 0.5d0 * (a43 * b34 - a33 - a44)
    r = t * t + a34 * a43 - a33 * a44
    if (r .lt. 0.0d0) goto 150
!   .......... determine single shift zeroth column of a ..........
    ish = 1
    r = dsqrt(r)
    sh = -t + r
    s = -t - r
    if (dabs(s-a44) .lt. dabs(sh-a44)) sh = s
!   .......... look for two consecutive small
!              sub-diagonal elements of a.
!              for l=en-2 step -1 until ld do -- ..........
    do ll = ld, enm2
        l = enm2 + ld - ll
        if (l .eq. ld) goto 140
        lm1 = l - 1
        l1 = l + 1
        t = a(l,l)
        if (dabs(b(l,l)) .gt. epsb) t = t - sh * b(l,l)
        if (dabs(a(l,lm1)) .le. dabs(t/a(l1,l)) * epsa) goto 100
    end do

140 a1 = a11 - sh
    a2 = a21
    if (l .ne. ld) a(l,lm1) = -a(l,lm1)
    goto 160
!   .......... determine double shift zeroth column of a ..........
150 a12 = a(l,l1) / b22
    a22 = a(l1,l1) / b22
    b12 = b(l,l1) / b22
    a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11) /&
         a21 + a12 - a11 * b12
    a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11) + a43 * b34
    a3 = a(l1+1,l1) / b22
    goto 160
!   .......... ad hoc shift ..........
155 a1 = 0.0d0
    a2 = 1.0d0
    a3 = 1.1605d0
160 its = its + 1
    itn = itn - 1
    if (.not. matz) lor1 = ld
!   .......... main loop ..........
    do k = l, na
        notlas = k .ne. na .and. ish .eq. 2
        k1 = k + 1
        k2 = k + 2
        km1 = max0(k-1,l)
        ll = min0(en,k1+ish)
        if (notlas) goto 190
!   .......... zero a(k+1,k-1) ..........
        if (k .eq. l) goto 170
        a1 = a(k,km1)
        a2 = a(k1,km1)
170     s = dabs(a1) + dabs(a2)
        if (s .eq. 0.0d0) goto 70
        u1 = a1 / s
        u2 = a2 / s
        r = dsign(dsqrt(u1*u1+u2*u2),u1)
        v1 = -(u1 + r) / r
        v2 = -u2 / r
        u2 = v2 / v1
        do j = km1, enorn
            t = a(k,j) + u2 * a(k1,j)
            a(k,j) = a(k,j) + t * v1
            a(k1,j) = a(k1,j) + t * v2
            t = b(k,j) + u2 * b(k1,j)
            b(k,j) = b(k,j) + t * v1
            b(k1,j) = b(k1,j) + t * v2
        end do
        if (k .ne. l) a(k1,km1) = 0.0d0
        goto 240
!   .......... zero a(k+1,k-1) and a(k+2,k-1) ..........
190     if (k .eq. l) goto 200
        a1 = a(k,km1)
        a2 = a(k1,km1)
        a3 = a(k2,km1)
200     s = dabs(a1) + dabs(a2) + dabs(a3)
        if (s .eq. 0.0d0) cycle
        u1 = a1 / s
        u2 = a2 / s
        u3 = a3 / s
        r = dsign(dsqrt(u1*u1+u2*u2+u3*u3),u1)
        v1 = -(u1 + r) / r
        v2 = -u2 / r
        v3 = -u3 / r
        u2 = v2 / v1
        u3 = v3 / v1
        do j = km1, enorn
            t = a(k,j) + u2 * a(k1,j) + u3 * a(k2,j)
            a(k,j) = a(k,j) + t * v1
            a(k1,j) = a(k1,j) + t * v2
            a(k2,j) = a(k2,j) + t * v3
            t = b(k,j) + u2 * b(k1,j) + u3 * b(k2,j)
            b(k,j) = b(k,j) + t * v1
            b(k1,j) = b(k1,j) + t * v2
            b(k2,j) = b(k2,j) + t * v3
        end do
        if (k .eq. l) goto 220
        a(k1,km1) = 0.0d0
        a(k2,km1) = 0.0d0
!   .......... zero b(k+2,k+1) and b(k+2,k) ..........
220     s = dabs(b(k2,k2)) + dabs(b(k2,k1)) + dabs(b(k2,k))
        if (s .eq. 0.0d0) goto 240
        u1 = b(k2,k2) / s
        u2 = b(k2,k1) / s
        u3 = b(k2,k) / s
        r = dsign(dsqrt(u1*u1+u2*u2+u3*u3),u1)
        v1 = -(u1 + r) / r
        v2 = -u2 / r
        v3 = -u3 / r
        u2 = v2 / v1
        u3 = v3 / v1
        do i = lor1, ll
            t = a(i,k2) + u2 * a(i,k1) + u3 * a(i,k)
            a(i,k2) = a(i,k2) + t * v1
            a(i,k1) = a(i,k1) + t * v2
            a(i,k) = a(i,k) + t * v3
            t = b(i,k2) + u2 * b(i,k1) + u3 * b(i,k)
            b(i,k2) = b(i,k2) + t * v1
            b(i,k1) = b(i,k1) + t * v2
            b(i,k) = b(i,k) + t * v3
        end do
        b(k2,k) = 0.0d0
        b(k2,k1) = 0.0d0
        if (.not. matz) goto 240
        do i = 1, n
            t = z(i,k2) + u2 * z(i,k1) + u3 * z(i,k)
            z(i,k2) = z(i,k2) + t * v1
            z(i,k1) = z(i,k1) + t * v2
            z(i,k) = z(i,k) + t * v3
        end do
!   .......... zero b(k+1,k) ..........
240     s = dabs(b(k1,k1)) + dabs(b(k1,k))
        if (s .eq. 0.0d0) cycle
        u1 = b(k1,k1) / s
        u2 = b(k1,k) / s
        r = dsign(dsqrt(u1*u1+u2*u2),u1)
        v1 = -(u1 + r) / r
        v2 = -u2 / r
        u2 = v2 / v1
        do i = lor1, ll
            t = a(i,k1) + u2 * a(i,k)
            a(i,k1) = a(i,k1) + t * v1
            a(i,k) = a(i,k) + t * v2
            t = b(i,k1) + u2 * b(i,k)
            b(i,k1) = b(i,k1) + t * v1
            b(i,k) = b(i,k) + t * v2
        end do
        b(k1,k) = 0.0d0
        if (.not. matz) cycle
        do i = 1, n
            t = z(i,k1) + u2 * z(i,k)
            z(i,k1) = z(i,k1) + t * v1
            z(i,k) = z(i,k) + t * v2
        end do
    end do
!   .......... end qz step ..........
    goto 70
!   .......... set error -- all eigenvalues have not
!              converged after 30*n iterations ..........
998 ierr = en
!     .......... save epsb for use by qzval and qzvec ..........
999 if (n .gt. 1) b(n,1) = epsb
end subroutine


subroutine qzval(nm,n,a,b,alfr,alfi,beta,matz,z)
    implicit none

    integer i,j,n,en,na,nm,nn,isw
    real(kind(1.0d0)) a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)
    real(kind(1.0d0)) c,d,e,r,s,t,an,a1,a2,bn,cq,cz,di,dr,ei,ti,tr,u1
    real(kind(1.0d0)) u2,v1,v2,a1i,a11,a12,a2i,a21,a22,b11,b12,b22,sqi,sqr
    real(kind(1.0d0)) ssi,ssr,szi,szr,a11i,a11r,a12i,a12r,a22i,a22r,epsb
    logical matz

!   this subroutine is the third step of the qz algorithm
!   for solving generalized matrix eigenvalue problems,
!   siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!   this subroutine accepts a pair of real matrices, one of them
!   in quasi-triangular form and the other in upper triangular form.
!   it reduces the quasi-triangular matrix further, so that any
!   remaining 2-by-2 blocks correspond to pairs of complex
!   eigenvalues, and returns quantities whose ratios give the
!   generalized eigenvalues.  it is usually preceded by  qzhes
!   and  qzit  and may be followed by  qzvec.
!
!   on input
!
!     nm must be set to the row dimension of two-dimensional
!     array parameters as declared in the calling program
!     dimension statement.
!
!     n is the order of the matrices.
!
!     a contains a real upper quasi-triangular matrix.
!
!     b contains a real upper triangular matrix.  in addition,
!     location b(n,1) contains the tolerance quantity (epsb)
!     computed and saved in  qzit.
!
!     matz should be set to .true. if the right hand transformations
!     are to be accumulated for later use in computing
!     eigenvectors, and to .false. otherwise.
!
!     z contains, if matz has been set to .true., the
!     transformation matrix produced in the reductions by qzhes
!     and qzit, if performed, or else the identity matrix.
!     if matz has been set to .false., z is not referenced.
!
!   on output
!
!     a has been reduced further to a quasi-triangular matrix
!     in which all nonzero subdiagonal elements correspond to
!     pairs of complex eigenvalues.
!
!     b is still in upper triangular form, although its elements
!     have been altered.  b(n,1) is unaltered.
!
!     alfr and alfi contain the real and imaginary parts of the
!     diagonal elements of the triangular matrix that would be
!     obtained if a were reduced completely to triangular form
!     by unitary transformations.  non-zero values of alfi occur
!     in pairs, the first member positive and the second negative.
!
!     beta contains the diagonal elements of the corresponding b,
!     normalized to be real and non-negative.  the generalized
!     eigenvalues are then the ratios ((alfr+i*alfi)/beta).
!
!     z contains the product of the right hand transformations
!     (for all three steps) if matz has been set to .true.
!
!   questions and comments should be directed to burton s. garbow,
!   mathematics and computer science div, argonne national laboratory
!
!   this version dated august 1983.
!
!   ------------------------------------------------------------------

    epsb = b(n,1)
    isw = 1
!   .......... find eigenvalues of quasi-triangular matrices.
!              for en=n step -1 until 1 do -- ..........
    do nn = 1, n
        en = n + 1 - nn
        na = en - 1
        if (isw .eq. 2) goto 505
        if (en .eq. 1) goto 410
        if (a(en,na) .ne. 0.0d0) goto 420
!   .......... 1-by-1 block, one real root ..........
410     alfr(en) = a(en,en)
        if (b(en,en) .lt. 0.0d0) alfr(en) = -alfr(en)
        beta(en) = dabs(b(en,en))
        alfi(en) = 0.0d0
        cycle !goto 510
!   .......... 2-by-2 block ..........
420     if (dabs(b(na,na)) .le. epsb) goto 455
        if (dabs(b(en,en)) .gt. epsb) goto 430
        a1 = a(en,en)
        a2 = a(en,na)
        bn = 0.0d0
        goto 435
430     an = dabs(a(na,na)) + dabs(a(na,en)) + dabs(a(en,na)) + dabs(a(en,en))
        bn = dabs(b(na,na)) + dabs(b(na,en)) + dabs(b(en,en))
        a11 = a(na,na) / an
        a12 = a(na,en) / an
        a21 = a(en,na) / an
        a22 = a(en,en) / an
        b11 = b(na,na) / bn
        b12 = b(na,en) / bn
        b22 = b(en,en) / bn
        e = a11 / b11
        ei = a22 / b22
        s = a21 / (b11 * b22)
        t = (a22 - e * b22) / b22
        if (dabs(e) .le. dabs(ei)) goto 431
        e = ei
        t = (a11 - e * b11) / b11
431     c = 0.5d0 * (t - s * b12)
        d = c * c + s * (a12 - e * b12)
        if (d .lt. 0.0d0) goto 480
!   .......... two real roots.
!              zero both a(en,na) and b(en,na) ..........
        e = e + (c + dsign(dsqrt(d),c))
        a11 = a11 - e * b11
        a12 = a12 - e * b12
        a22 = a22 - e * b22
        if (dabs(a11) + dabs(a12) .lt. dabs(a21) + dabs(a22)) goto 432
        a1 = a12
        a2 = a11
        goto 435
432     a1 = a22
        a2 = a21
!   .......... choose and apply real z ..........
435     s = dabs(a1) + dabs(a2)
        u1 = a1 / s
        u2 = a2 / s
        r = dsign(dsqrt(u1*u1+u2*u2),u1)
        v1 = -(u1 + r) / r
        v2 = -u2 / r
        u2 = v2 / v1
        do i = 1, en
            t = a(i,en) + u2 * a(i,na)
            a(i,en) = a(i,en) + t * v1
            a(i,na) = a(i,na) + t * v2
            t = b(i,en) + u2 * b(i,na)
            b(i,en) = b(i,en) + t * v1
            b(i,na) = b(i,na) + t * v2
        end do

        if (.not. matz) goto 450
        do i = 1, n
            t = z(i,en) + u2 * z(i,na)
            z(i,en) = z(i,en) + t * v1
            z(i,na) = z(i,na) + t * v2
        end do
450     if (bn .eq. 0.0d0) goto 475
        if (an .lt. dabs(e) * bn) goto 455
        a1 = b(na,na)
        a2 = b(en,na)
        goto 460
455     a1 = a(na,na)
        a2 = a(en,na)
!   .......... choose and apply real q ..........
460     s = dabs(a1) + dabs(a2)
        if (s .eq. 0.0d0) goto 475
        u1 = a1 / s
        u2 = a2 / s
        r = dsign(dsqrt(u1*u1+u2*u2),u1)
        v1 = -(u1 + r) / r
        v2 = -u2 / r
        u2 = v2 / v1
        do j = na, n
            t = a(na,j) + u2 * a(en,j)
            a(na,j) = a(na,j) + t * v1
            a(en,j) = a(en,j) + t * v2
            t = b(na,j) + u2 * b(en,j)
            b(na,j) = b(na,j) + t * v1
            b(en,j) = b(en,j) + t * v2
        end do
475     a(en,na) = 0.0d0
        b(en,na) = 0.0d0
        alfr(na) = a(na,na)
        alfr(en) = a(en,en)
        if (b(na,na) .lt. 0.0d0) alfr(na) = -alfr(na)
        if (b(en,en) .lt. 0.0d0) alfr(en) = -alfr(en)
        beta(na) = dabs(b(na,na))
        beta(en) = dabs(b(en,en))
        alfi(en) = 0.0d0
        alfi(na) = 0.0d0
        goto 505
!   .......... two complex roots ..........
480     e = e + c
        ei = dsqrt(-d)
        a11r = a11 - e * b11
        a11i = ei * b11
        a12r = a12 - e * b12
        a12i = ei * b12
        a22r = a22 - e * b22
        a22i = ei * b22
        if (dabs(a11r) + dabs(a11i) + dabs(a12r) + dabs(a12i) .lt.&
           dabs(a21) + dabs(a22r) + dabs(a22i)) goto 482
        a1 = a12r
        a1i = a12i
        a2 = -a11r
        a2i = -a11i
        goto 485
482     a1 = a22r
        a1i = a22i
        a2 = -a21
        a2i = 0.0d0
!   .......... choose complex z ..........
485     cz = dsqrt(a1*a1+a1i*a1i)
        if (cz .eq. 0.0d0) goto 487
        szr = (a1 * a2 + a1i * a2i) / cz
        szi = (a1 * a2i - a1i * a2) / cz
        r = dsqrt(cz*cz+szr*szr+szi*szi)
        cz = cz / r
        szr = szr / r
        szi = szi / r
        goto 490
487     szr = 1.0d0
        szi = 0.0d0
490     if (an .lt. (dabs(e) + ei) * bn) goto 492
        a1 = cz * b11 + szr * b12
        a1i = szi * b12
        a2 = szr * b22
        a2i = szi * b22
        goto 495
492     a1 = cz * a11 + szr * a12
        a1i = szi * a12
        a2 = cz * a21 + szr * a22
        a2i = szi * a22
!   .......... choose complex q ..........
495     cq = dsqrt(a1*a1+a1i*a1i)
        if (cq .eq. 0.0d0) goto 497
        sqr = (a1 * a2 + a1i * a2i) / cq
        sqi = (a1 * a2i - a1i * a2) / cq
        r = dsqrt(cq*cq+sqr*sqr+sqi*sqi)
        cq = cq / r
        sqr = sqr / r
        sqi = sqi / r
        goto 500
497     sqr = 1.0d0
        sqi = 0.0d0
!   .......... compute diagonal elements that would result
!              if transformations were applied ..........
500     ssr = sqr * szr + sqi * szi
        ssi = sqr * szi - sqi * szr
        i = 1
        tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22
        ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22
        dr = cq * cz * b11 + cq * szr * b12 + ssr * b22
        di = cq * szi * b12 + ssi * b22
        goto 503
502     i = 2
        tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22
        ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21
        dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22
        di = -ssi * b11 - sqi * cz * b12
503     t = ti * dr - tr * di
        j = na
        if (t .lt. 0.0d0) j = en
        r = dsqrt(dr*dr+di*di)
        beta(j) = bn * r
        alfr(j) = an * (tr * dr + ti * di) / r
        alfi(j) = an * t / r
        if (i .eq. 1) goto 502
505     isw = 3 - isw
    end do
    b(n,1) = epsb
end subroutine


subroutine qzvec(nm,n,a,b,alfr,alfi,beta,z)
    implicit none

    integer i,j,k,m,n,en,ii,jj,na,nm,nn,isw,enm2
    real(kind(1.0d0)) a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)
    real(kind(1.0d0)) d,q,r,s,t,w,x,y,di,dr,ra,rr,sa,ti,tr,t1,t2,w1,x1
    real(kind(1.0d0)) zz,z1,alfm,almi,almr,betm,epsb

!   this subroutine is the optional fourth step of the qz algorithm
!   for solving generalized matrix eigenvalue problems,
!   siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!   this subroutine accepts a pair of real matrices, one of them in
!   quasi-triangular form (in which each 2-by-2 block corresponds to
!   a pair of complex eigenvalues) and the other in upper triangular
!   form.  it computes the eigenvectors of the triangular problem and
!   transforms the results back to the original coordinate system.
!   it is usually preceded by  qzhes,  qzit, and  qzval.
!
!   on input
!
!     nm must be set to the row dimension of two-dimensional
!     array parameters as declared in the calling program
!     dimension statement.
!
!     n is the order of the matrices.
!
!     a contains a real upper quasi-triangular matrix.
!
!     b contains a real upper triangular matrix.  in addition,
!     location b(n,1) contains the tolerance quantity (epsb)
!     computed and saved in  qzit.
!
!     alfr, alfi, and beta  are vectors with components whose
!     ratios ((alfr+i*alfi)/beta) are the generalized
!     eigenvalues.  they are usually obtained from  qzval.
!
!     z contains the transformation matrix produced in the
!     reductions by  qzhes,  qzit, and  qzval, if performed.
!     if the eigenvectors of the triangular problem are
!     desired, z must contain the identity matrix.
!
!   on output
!
!     a is unaltered.  its subdiagonal elements provide information
!     about the storage of the complex eigenvectors.
!
!     b has been destroyed.
!
!     alfr, alfi, and beta are unaltered.
!
!     z contains the real and imaginary parts of the eigenvectors.
!     if alfi(i) .eq. 0.0, the i-th eigenvalue is real and
!     the i-th column of z contains its eigenvector.
!     if alfi(i) .ne. 0.0, the i-th eigenvalue is complex.
!     if alfi(i) .gt. 0.0, the eigenvalue is the first of
!     a complex pair and the i-th and (i+1)-th columns
!     of z contain its eigenvector.
!     if alfi(i) .lt. 0.0, the eigenvalue is the second of
!     a complex pair and the (i-1)-th and i-th columns
!     of z contain the conjugate of its eigenvector.
!     each eigenvector is normalized so that the modulus
!     of its largest component is 1.0 .
!
!   questions and comments should be directed to burton s. garbow,
!   mathematics and computer science div, argonne national laboratory
!
!   this version dated august 1983.
!
!   ------------------------------------------------------------------

    epsb = b(n,1)
    isw = 1
!   .......... for en=n step -1 until 1 do -- ..........
    do nn = 1, n
        en = n + 1 - nn
        na = en - 1
        if (isw .eq. 2) goto 795
        if (alfi(en) .ne. 0.0d0) goto 710
!   .......... real vector ..........
        m = en
        b(en,en) = 1.0d0
        if (na .eq. 0) cycle
        alfm = alfr(m)
        betm = beta(m)
!   .......... for i=en-1 step -1 until 1 do -- ..........
        do ii = 1, na
            i = en - ii
            w = betm * a(i,i) - alfm * b(i,i)
            r = 0.0d0
            do j = m, en
                r = r + (betm * a(i,j) - alfm * b(i,j)) * b(j,en)
            end do
            if (i .eq. 1 .or. isw .eq. 2) goto 630
            if (betm * a(i,i-1) .eq. 0.0d0) goto 630
            zz = w
            s = r
            goto 690
630         m = i
            if (isw .eq. 2) goto 640
!   .......... real 1-by-1 block ..........
            t = w
            if (w .eq. 0.0d0) t = epsb
            b(i,en) = -r / t
            cycle
!   .......... real 2-by-2 block ..........
640         x = betm * a(i,i+1) - alfm * b(i,i+1)
            y = betm * a(i+1,i)
            q = w * zz - x * y
            t = (x * s - zz * r) / q
            b(i,en) = t
            if (dabs(x) .le. dabs(zz)) goto 650
            b(i+1,en) = (-r - w * t) / x
            goto 690
650         b(i+1,en) = (-s - y * t) / zz
690         isw = 3 - isw
        end do
!     .......... end real vector ..........
        cycle
!     .......... complex vector ..........
710     m = na
        almr = alfr(m)
        almi = alfi(m)
        betm = beta(m)
!   .......... last vector component chosen imaginary so that
!              eigenvector matrix is triangular ..........
        y = betm * a(en,na)
        b(na,na) = -almi * b(en,en) / y
        b(na,en) = (almr * b(en,en) - betm * a(en,en)) / y
        b(en,na) = 0.0d0
        b(en,en) = 1.0d0
        enm2 = na - 1
        if (enm2 .eq. 0) goto 795
!   .......... for i=en-2 step -1 until 1 do -- ..........
        do ii = 1, enm2
            i = na - ii
            w = betm * a(i,i) - almr * b(i,i)
            w1 = -almi * b(i,i)
            ra = 0.0d0
            sa = 0.0d0
            do j = m, en
                x = betm * a(i,j) - almr * b(i,j)
                x1 = -almi * b(i,j)
                ra = ra + x * b(j,na) - x1 * b(j,en)
                sa = sa + x * b(j,en) + x1 * b(j,na)
            end do
            if (i .eq. 1 .or. isw .eq. 2) goto 770
            if (betm * a(i,i-1) .eq. 0.0d0) goto 770
            zz = w
            z1 = w1
            r = ra
            s = sa
            isw = 2
            cycle
770         m = i
            if (isw .eq. 2) goto 780
!   .......... complex 1-by-1 block ..........
            tr = -ra
            ti = -sa
773         dr = w
            di = w1
!   .......... complex divide (t1,t2) = (tr,ti) / (dr,di) ..........
775         if (dabs(di) .gt. dabs(dr)) goto 777
            rr = di / dr
            d = dr + di * rr
            t1 = (tr + ti * rr) / d
            t2 = (ti - tr * rr) / d
            if (isw .ne. 2) then
                goto 787
            else
                goto 782
            end if
777         rr = dr / di
            d = dr * rr + di
            t1 = (tr * rr + ti) / d
            t2 = (ti * rr - tr) / d
            if (isw .ne. 2) then
                goto 787
            else
                goto 782
            end if
!   .......... complex 2-by-2 block ..........
780         x = betm * a(i,i+1) - almr * b(i,i+1)
            x1 = -almi * b(i,i+1)
            y = betm * a(i+1,i)
            tr = y * ra - w * r + w1 * s
            ti = y * sa - w * s - w1 * r
            dr = w * zz - w1 * z1 - x * y
            di = w * z1 + w1 * zz - x1 * y
            if (dr .eq. 0.0d0 .and. di .eq. 0.0d0) dr = epsb
            goto 775
782         b(i+1,na) = t1
            b(i+1,en) = t2
            isw = 1
            if (dabs(y) .gt. dabs(w) + dabs(w1)) goto 785
            tr = -ra - x * b(i+1,na) + x1 * b(i+1,en)
            ti = -sa - x * b(i+1,en) - x1 * b(i+1,na)
            goto 773
785         t1 = (-r - zz * b(i+1,na) + z1 * b(i+1,en)) / y
            t2 = (-s - zz * b(i+1,en) - z1 * b(i+1,na)) / y
787         b(i,na) = t1
            b(i,en) = t2
        end do
!   .......... end complex vector ..........
795     isw = 3 - isw
    end do
!   .......... end back substitution.
!              transform to original coordinate system.
!              for j=n step -1 until 1 do -- ..........
    do jj = 1, n
        j = n + 1 - jj
        do i = 1, n
            zz = 0.0d0
            do k = 1, j
                zz = zz + z(i,k) * b(k,j)
            end do 
            z(i,j) = zz
        end do
    end do
!   .......... normalize so that modulus of largest
!              component of each vector is 1.
!              (isw is 1 initially from before) ..........
    do j = 1, n
        d = 0.0d0
        if (isw .eq. 2) goto 920
        if (alfi(j) .ne. 0.0d0) goto 945
        do i = 1, n
            if (dabs(z(i,j)) .gt. d) d = dabs(z(i,j))
        end do
        do i = 1, n
            z(i,j) = z(i,j) / d
        end do
        cycle
920     do i = 1, n
            r = dabs(z(i,j-1)) + dabs(z(i,j))
            if (r .ne. 0.0d0) r = r * dsqrt((z(i,j-1)/r)**2 + (z(i,j)/r)**2)
            if (r .gt. d) d = r
        end do
        do i = 1, n
            z(i,j-1) = z(i,j-1) / d
            z(i,j) = z(i,j) / d
        end do
945     isw = 3 - isw
    end do
end subroutine


subroutine rgg(nm,n,a,b,alfr,alfi,beta,matz,z,ierr)
    implicit none

    integer n,nm,ierr,matz
    real(kind(1.0d0)) a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)
    logical tf

!   this subroutine calls the recommended sequence of
!   subroutines from the eigensystem subroutine package (eispack)
!   to find the eigenvalues and eigenvectors (if desired)
!   for the real general generalized eigenproblem  ax = (lambda)bx.
!
!   on input
!
!     nm  must be set to the row dimension of the two-dimensional
!     array parameters as declared in the calling program
!     dimension statement.
!
!     n  is the order of the matrices  a  and  b.
!
!     a  contains a real general matrix.
!
!     b  contains a real general matrix.
!
!     matz  is an integer variable set equal to zero if
!     only eigenvalues are desired.  otherwise it is set to
!     any non-zero integer for both eigenvalues and eigenvectors.
!
!   on output
!
!     alfr  and  alfi  contain the real and imaginary parts,
!     respectively, of the numerators of the eigenvalues.
!
!     beta  contains the denominators of the eigenvalues,
!     which are thus given by the ratios  (alfr+i*alfi)/beta.
!     complex conjugate pairs of eigenvalues appear consecutively
!     with the eigenvalue having the positive imaginary part first.
!
!     z  contains the real and imaginary parts of the eigenvectors
!     if matz is not zero.  if the j-th eigenvalue is real, the
!     j-th column of  z  contains its eigenvector.  if the j-th
!     eigenvalue is complex with positive imaginary part, the
!     j-th and (j+1)-th columns of  z  contain the real and
!     imaginary parts of its eigenvector.  the conjugate of this
!     vector is the eigenvector for the conjugate eigenvalue.
!
!     ierr  is an integer output variable set equal to an error
!     completion code described in the documentation for qzit.
!     the normal completion code is zero.
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
    else
        if (matz .eq. 0) then
!   .......... find eigenvalues only ..........
            tf = .false.
            call  qzhes(nm,n,a,b,tf,z)
            call  qzit(nm,n,a,b,0.0d0,tf,z,ierr)
            call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z)
            return
        else
!   .......... find both eigenvalues and eigenvectors ..........
            tf = .true.
            call  qzhes(nm,n,a,b,tf,z)
            call  qzit(nm,n,a,b,0.0d0,tf,z,ierr)
            call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z)
            if (ierr .ne. 0) return
            call  qzvec(nm,n,a,b,alfr,alfi,beta,z)
        endif
    endif
end subroutine
