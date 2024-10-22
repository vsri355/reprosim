module solve
!*Description:* This module contains solvers required for placenta problems
!
! Solvers included in this module are:
!                                                                          
! BICGSTAB == The STABalized BIConjugate Gradient method
!                                                                           
! GMRES == Generalised Minimal RESidual method
!
! GMRES solver. The following is the code from mgmres.f90, downloaded from
! www.people.sc.fsu.edu/~jbukhardt/f_src/mgmres/mgmres.html on 13/01/2016
!
! Licensing:
!    This code is distributed under the GNU LGPL license.
! Modified:
!    17 July 2007
! Author:
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
! Reference:
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!    Input, integer :: N, the order of the system.
!
!    Input, integer :: NZ_NUM, the number of nonzeros.
!
!    Input, integer :: IA(N+1), JA(NZ_NUM), the row and column
!    indices of the matrix values.  The row vector has been compressed.
!
!    Input, real(dp) :: A(NZ_NUM), the matrix values.
!
!    Input, real(dp) :: X(N), the vector to be multiplied by A'.
!
!    Output, real(dp) :: W(N), the value of A'*X.
!
! Information:
!    The Sparse Compressed Row storage format is used.
!
!    The matrix A is assumed to be sparse.  To save on storage, only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA[I] through IA[I+1]-1.
!
!       
  use math_utilities
  implicit none
  private
  public BICGSTAB_LinSolv
  public pmgmres_ilu_cr
contains
!
!#######################################################################
!

subroutine BICGSTAB_LinSolv(MatrixSize,NonZeros,RHS,Solution,SparseCol,SparseRow,SparseVal,TOLER,MaxIter)
!*Description:*
!   This subroutine uses a preconditioned biconjugate gradient method to solve a
!   linear system Ax=b. Currently this solver is hard-coded to use the diaganol
!   entries of A as a preconditioner (the Jacobi preconditioner). However, this
!   could easily be added as an input variable. Sparsity is defined by compressed
!   row storage.
! Created by ARC: 15-12-2015
!
!   Input:
!   MatrixSize -> the size of the matrix.
!   NonZeros -> The number of non-zero entries.
!   SparseCol -> Define the column indices that have non-zero entries.
!   SparseRow ->Define the indices of SparseCol that represent each row.
!   SparseVal -> The values of non-zero entries.
!   RHS -> The vector on the RHS of the linear system
!
!   FLAGS - 0=Solution found to tolerance
!           1=No convergence max iterations reached
!           -1=breakdown rho=0
!           -2=breakdown omega=0

    use arrays

    !Input/Output Variables
    integer, intent(in) :: MatrixSize,NonZeros
    integer, intent(in) :: SparseCol(NonZeros),SparseRow(MatrixSize+1),MaxIter
    double precision :: Solution(MatrixSize),SparseVal(NonZeros),RHS(MatrixSize),TOLER

    !Local Variables
    integer :: i,k,j,kstart,kend,Iter,FLAG
    double precision :: Ax(MatrixSize),bnrm2,error_r,PreCon(MatrixSize),PreConInv(MatrixSize)
    double precision :: r(MatrixSize),rtilde(MatrixSize),alpha,beta,omega,rho,rho_1
    double precision :: p(MatrixSize),v(MatrixSize),phat(MatrixSize),s(MatrixSize),snorm,resid,shat(MatrixSize),t(MatrixSize)
    real :: start,finish
    call CPU_TIME(start)
    print *, 'In Solver, Tolerance=', TOLER, 'Max iterations=', MaxIter

! PRECONDITIONER MATRIX - THE DIAGONALS OF THE MATRIX A P=diag(A), Jacobi preconditioner Inverse is deltaij/Aij
    DO i=1,MatrixSize
        kstart=SparseRow(i)
        kend=SparseRow(i+1)-1
        DO k=kstart,kend
            IF(SparseCol(k).eq.i)THEN
                PreCon(i)=SparseVal(k)
                PreConInv(i)=1/PreCon(i)
            ENDIF
        ENDDO
    ENDDO


    !bnrm2=norm(RHS)
    bnrm2=sqrt(DOT_PRODUCT(RHS,RHS))!Checked against matlab
    IF(bnrm2.eq.0.0_dp) bnrm2=1.0_dp


    error_r=0.0_dp
    !CALL ax_cr(MatrixSize,SparseRow,SparseCol,SparseVal,solution,Ax)
    DO i=1,MatrixSize !This loop depends on sparsity structure
        Ax(i)=0.0_dp !Initialise Ax for this entry
        kstart=SparseRow(i)
        kend=SparseRow(i+1)-1
        DO k=kstart,kend
            Ax(i)=Ax(i)+SparseVal(k)*solution(SparseCol(k)) !A*X
        ENDDO
    ENDDO
    r=RHS-Ax
    rtilde=r

    error_r=sqrt(DOT_PRODUCT(r,r))/bnrm2 !Checked against matlab


    IF(error_r.lt.TOLER)THEN
        write(*,*) 'Initial guess is perfect, Exiting solver',TOLER
        FLAG=0
    RETURN
    ENDIF
    alpha=0.0_dp
    beta=0.0_dp
    omega=1.0_dp

!START OF ITERATIVE LOOP
     iterative_loop: DO Iter=1,MaxIter
        !rho is dot product of r and t tilde
        rho=DOT_PRODUCT(r,rtilde) !checked against matlab
        IF(rho.EQ.0.0_dp) EXIT iterative_loop
        IF(ITER.GT.1.0_dp)THEN
            beta=(rho/rho_1)*(alpha/omega)
            p=r+beta*(p-omega*v)
            phat=PreConInv*p
            v=0.0_dp*v
        ELSE
            p=r
            phat=PreConInv*p
            v=0.0_dp*p
        ENDIF


        !v=ApHat
        DO i=1,MatrixSize !This loop depends on sparsity structure
            kstart=SparseRow(i)
            kend=SparseRow(i+1)-1
            DO k=kstart,kend
                v(i)=v(i)+SparseVal(k)*phat(SparseCol(k))
            ENDDO
        ENDDO
        !v is correct against matlab

        alpha=rho/DOT_PRODUCT(rtilde,v) !alpha is correct first iter
        s=r-alpha*v

        snorm=sqrt(DOT_PRODUCT(s,s))!snorm correct against matlab

        !!EARLY CONVERGENCE TEST
        IF(snorm.LE.TOLER)THEN
            solution=solution+alpha*phat
            resid=snorm/bnrm2
            EXIT iterative_loop
        ENDIF


        !Stabiliser shat=precon-1*s
        shat=PreConInv*s
        t=0.0_dp*s
        DO i=1,MatrixSize
            kstart=SparseRow(i)
            kend=SparseRow(i+1)-1
            DO k=kstart,kend
                t(i)=t(i)+SparseVal(k)*shat(SparseCol(k)) !A*X !correct first iter
            ENDDO
        ENDDO
        omega=DOT_PRODUCT(t,s)/DOT_PRODUCT(t,t)!omega corect first iter


        solution=solution+alpha*phat+omega*shat
        r=s-omega*t
        error_r=SQRT(DOT_PRODUCT(r,r))/bnrm2

        IF(error_r.LE.TOLER) EXIT iterative_loop
        IF(omega.EQ.0.0_dp) EXIT iterative_loop
        rho_1=rho!

    ENDDO iterative_loop
    IF(error_r.LE.TOLER.OR.snorm.LE.TOLER)THEN
        IF(snorm.LE.TOLER)THEN
            error_r=SQRT(DOT_PRODUCT(s,s))/bnrm2
        ENDIF
!        FLAG=0
        print *, 'Converged with error=',error_r
    ELSEIF(omega.EQ.0)THEN
!        FLAG=-2
        print *, 'Not converged omega =0'
    ELSEIF(rho.EQ.0)THEN
!        FLAG=-1
        print *, 'Not converged rho =0'
    ELSE
!        FLAG=1
        print *, 'Not converged max no of iterations reached'
        print *, 'Error=',error_r
    ENDIF
    CALL CPU_TIME(finish)
    print '("Time= ",f12.9,"seconds")',finish-start



end subroutine BICGSTAB_LinSolv
!
!#############################################################################
!
subroutine pmgmres_ilu_cr (n,nz_num,ia,ja,a, x, rhs, itr_max, mr, tol_abs, tol_rel,FLAG)
! PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
!    This routine uses the incomplete LU decomposition for the
!    preconditioning.  This preconditioner requires that the sparse
!    matrix data structure supplies a storage position for each diagonal
!    element of the matrix A, and that each diagonal element of the
!    matrix A is not zero.
!    Input/output, real(dp) :: X(N); on input, an approximation to
!    the solution.  On output, an improved approximation.
!
!    Input, real(dp) :: RHS(N), the right hand side of the linear system.
!
!    Input, integer :: ITR_MAX, the maximum number of (outer)
!    iterations to take.
!
!    Input, integer :: MR, the maximum number of (inner) iterations
!    to take.  MR must be less than N.
!
!    Input, real(dp) :: TOL_ABS, an absolute tolerance applied to the
!    current residual.
!
!    Input, real(dp) :: TOL_REL, a relative tolerance comparing the
!    current residual to the initial residual.
!
    use arrays
    use diagnostics, only: enter_exit,get_diagnostics_level
    implicit none

    integer, intent(in) :: n !n=the size of the matrix
    !integer :: mr
    integer, intent(in) :: nz_num !Number of non-zeros within the matrix
    integer, intent(in) :: ia(n+1) !SparseRow
    integer, intent(in) :: ja(nz_num) !SparseCol
    real(dp), intent(in) :: a(nz_num) !sparseval
    real(dp), intent(inout) :: x(n) !Solution
    real(dp), intent(in) :: rhs(n) !rhs(n)

    integer, intent(in) :: itr_max
    integer, intent(in) :: mr
    real(dp), intent(in) :: tol_abs
    real(dp), intent(in) :: tol_rel
    integer, intent(out) :: FLAG

    !! Local variable
    real(dp) :: av
    real(dp), allocatable :: c(:)!(mr+1)
    real(dp) :: delta = 1.0d-03

    real(dp) :: htmp
    integer :: i
    integer :: itr
    integer :: itr_used
    integer :: j
    integer :: k
    integer :: k_copy

    real(dp), allocatable :: l(:)!l(ia(n+1)+1)
    integer :: l_size
    real(dp) :: mu
    real(dp), allocatable :: r(:)
    real(dp) :: rho
    real(dp) :: rho_tol
    real(dp), allocatable :: s(:)! (mr+1)
    integer :: ua(n)
    real(dp), allocatable :: v(:,:)
    real(dp), allocatable :: g(:)!(mr+1)
    real(dp), allocatable :: h(:,:)!(mr+1,mr)
    real(dp), allocatable :: y(:)!(mr+1)

    integer :: diagnostics_level
    character(len=60) :: sub_name


    sub_name = 'pmgmres_ilu_cr'
    call enter_exit(sub_name,1)
    call get_diagnostics_level(diagnostics_level)

    !Zero the number of local iterations
    itr_used = 0
    !Solver flag
    FLAG=0 !not converged

    call rearrange_cr( n, nz_num, ia, ja, a )

    call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )
    !!! Allocating memory to local parameters
    if(allocated(l)) deallocate (l)
    l_size = ia(n+1)+1
    allocate (l(l_size))
    if(allocated(r)) deallocate(r)
    allocate(r(n))
    if(allocated(v)) deallocate(v)
    allocate(v(n,mr+1))
    if(allocated(g)) deallocate(g)
    allocate(g(mr+1))
    if(allocated(h)) deallocate(h)
    allocate(h(mr+1,mr))
    if(allocated(c)) deallocate(c)
    allocate(c(mr+1))
    if(allocated(s)) deallocate(s)
    allocate(s(mr+1))
    if(allocated(y)) deallocate(y)
    allocate(y(mr+1))

    call ilu_cr ( n, nz_num, l_size, ia, ja, a, ua, l)

    if (diagnostics_level.GT.1) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'PMGMRES_ILU_CR'
       write ( *, '(a,i8)' ) '  Number of unknowns = ', n
    end if

    do itr = 1, itr_max !DO ITR

       call ax_cr ( n, nz_num, ia, ja, a, x, r )

       r(1:n) = rhs(1:n) - r(1:n)

       call lus_cr ( n, nz_num, l_size, ia, ja, l, ua, r)

       rho = sqrt ( dot_product ( r, r ) )

       if (diagnostics_level.GT.1) then
          write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
       end if

       if ( itr == 1 ) then
          rho_tol = rho * tol_rel
       end if

       v(1:n,1) = r(1:n) / rho

       g(1) = rho
       g(2:mr+1) = 0.0_dp

       h(1:mr+1,1:mr) = 0.0_dp

       do k = 1, mr !DO K

          k_copy = k

          call ax_cr ( n,nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) )

          call lus_cr ( n,nz_num, l_size, ia, ja, l, ua, v(1:n,k+1))

          av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

          do j = 1, k
             h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
             v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
          end do

          h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

          if ( ( av + delta * h(k+1,k)) == av ) then
             do j = 1, k
                htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
                h(j,k) = h(j,k) + htmp
                v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
             end do
             h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
          end if

          if ( h(k+1,k).ne.0.0_dp ) then
             v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
          end if

          if ( 1 < k ) then
             y(1:k+1) = h(1:k+1,k)
             do j = 1, k - 1
                call mult_givens ( c(j), s(j), j, y )
             end do
             h(1:k+1,k) = y(1:k+1)
          end if

          mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

          c(k) = h(k,k) / mu
          s(k) = -h(k+1,k) / mu
          h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
          h(k+1,k) = 0.0_dp
          call mult_givens ( c(k), s(k), k, g )

          rho = abs ( g(k+1) )

          itr_used = itr_used + 1

          if (diagnostics_level.GT.1) then
             write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
          end if

          if ( rho <= rho_tol .and. rho <= tol_abs ) then
             FLAG=1
             exit
          elseif (isnan(rho))then
             FLAG=2
             exit
          end if

       end do !ENDDO K

      k = k_copy - 1

      y(k+1) = g(k+1) / h(k+1,k+1)

      do i = k, 1, -1 !DO K BACK
         y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
      end do !ENDDO K BACK

      do i = 1, n !DO N
         x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
      end do !ENDDO n

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
          FLAG=1
         exit
      elseif (isnan(rho))then
         FLAG=2
         exit
      end if

    end do !ENDDO ITR
    
    if (diagnostics_level.GT.1) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
       write ( *, '(a,i6)' ) '  Iterations = ', itr_used
       write ( *, '(a,g14.6)' ) '  Final residual = ', rho
    end if
    !!! Deallocate local memory
    deallocate(l)
    deallocate(r)
    deallocate(v)
    deallocate(g)
    deallocate(h)
    deallocate(c)
    deallocate(s)
    deallocate(y)
    call enter_exit(sub_name,2)

end subroutine pmgmres_ilu_cr



end module solve
