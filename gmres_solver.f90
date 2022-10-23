module gmres_solver_mod
  !
  use dtype_mod
  !
  implicit none
  !
  private
  public :: gmres_solver_t
  !
  type :: gmres_solver_t
    ! GMRES dimension
    integer :: m
    ! The dimension of solution
    integer :: n
    ! Number of non-zero in the sparse matrix A
    ! Possibly larger than the real number of non-zero in A
    integer :: nnz
    ! The solution
    real(wp), dimension(:), allocatable :: u
    ! RHS
    real(wp), dimension(:), allocatable :: b
    real(wp), dimension(:), allocatable :: u_exact
    !
    real(wp), dimension(:), allocatable :: A
    integer, dimension(:), allocatable :: ja
    integer, dimension(:), allocatable :: ia
    !
    real(wp), dimension(:, :), allocatable :: hess, hess_tmp
    real(wp), dimension(:, :), allocatable :: kryl
    !
    real(wp), dimension(:), allocatable :: work_w, work_z, work_y, dgelsy_work
    integer :: dgelsy_lwork
    integer, dimension(:), allocatable :: dgelsy_piv
    !
    integer :: gmres_pc_type
    integer :: rst_max
    real(wp) :: atol, rtol, dtol
    !> Instruction flag for GMRES. Specified by the user.
    !> 1: print the residual at each restart iteration
    integer :: gmres_info
    integer :: its
    !
    real(wp) :: idx2, idiag
    character(len=64) :: res_fname
    !
  contains
    !
    procedure :: SetExactSolution, SetB, AssembleA, PCSolve
    procedure :: Solve, ShowSaveResult, Free
    !
  end type gmres_solver_t
  !
  interface gmres_solver_t
    module procedure ConstructGMRESSolver
  end interface gmres_solver_t
  !
contains
  !
  function ConstructGMRESSolver(m, n, gmres_pc_type, rst_max, &
                                rtol, atol, dtol, gmres_info) result(slv)
    !>
    integer, intent(in) :: m
    !>
    integer, intent(in) :: n
    !>
    integer, intent(in) :: gmres_pc_type
    !>
    integer, intent(in) :: rst_max
    !>
    real(wp), intent(in) :: rtol
    !>
    real(wp), intent(in) :: atol
    !>
    real(wp), intent(in) :: dtol
    !>
    integer, intent(in) :: gmres_info
    !>
    type(gmres_solver_t) :: slv
    !
    integer :: i, nnz_per_row
    real(wp) :: dx
    !
  continue
    !
    slv%m = m
    slv%n = n
    slv%gmres_pc_type = gmres_pc_type
    slv%rst_max = rst_max
    slv%rtol = rtol
    slv%atol = atol
    slv%dtol = dtol
    slv%gmres_info = gmres_info
    !
    dx = one / real(slv%n-1, kind(dx))
    slv%idx2 = one / dx**2
    slv%idiag = -one / ( two * slv%idx2 )
    !
    nnz_per_row = 3
    slv%nnz = nnz_per_row * slv%n
    ! The solution is also initilized to zero here
    ! The initialization is required by GMRES
    allocate(slv%u(1:slv%n), source=zero)
    allocate(slv%b(1:slv%n), source=zero)
    allocate(slv%u_exact(1:slv%n), source=zero)
    allocate(slv%A(1:slv%nnz), source=zero)
    allocate(slv%ja(1:slv%nnz), source=0)
    allocate(slv%ia(1:slv%n+1), source=0)
    allocate(slv%work_w(1:slv%n), source=zero)
    allocate(slv%work_z(1:slv%n), source=zero)
    allocate(slv%work_y(1:slv%m+1), source=zero)
    ! For the value of dgelsy_lwork, check the API of LAPACK.
    slv%dgelsy_lwork = 4*m+1
    ! slv%dgelsy_lwork = 1024*m+1
    allocate(slv%dgelsy_work(1:slv%dgelsy_lwork), source=zero)
    allocate(slv%dgelsy_piv(1:slv%m), source=0)
    !
    allocate(slv%hess(1:slv%m+1, 1:slv%m), source=zero)
    allocate(slv%hess_tmp(1:slv%m+1, 1:slv%m), source=zero)
    allocate(slv%kryl(1:slv%n, 1:slv%m+1), source=zero)
    !
    call slv%SetB()
    call slv%SetExactSolution()
    !
    call slv%AssembleA()
    !
    select case (slv%gmres_pc_type)
    case (GMRES_No_PC)
      slv%res_fname = "gmres_no_pc.dat"
    case (GMRES_Jacobi)
      slv%res_fname = "gmres_jacobi.dat"
    case (GMRES_MG)
      slv%res_fname = "gmres_mg.dat"
    case default
      write(*, *) "gmres_pc_type invalid!"
      stop
    end select
    !
  end function ConstructGMRESSolver
  !
  subroutine SetB(this, r)
    !>
    class(gmres_solver_t) :: this
    !
    real(wp), dimension(:), intent(in), optional :: r
    !
    integer :: ierr, i
    real(wp) :: dx
    !
  continue
    !
    if (present(r)) then
      this%b(:) = r
    else
      dx = one / real(this%n-1, kind(dx))
      do i = 2, this%n-1
        this%b(i) = six * dx * real(i-1, kind(dx))
      end do
      this%b(1) = -one
      this%b(this%n) = zero
    end if
    !
  end subroutine SetB
  !
  subroutine SetExactSolution(this)
    !>
    class(gmres_solver_t) :: this
    !
    integer :: ierr, i
    real(wp) :: dx
    !
  continue
    dx = one / real(this%n-1, kind(dx))
    ! Analytic solution
    do i = 1, this%n
      this%u_exact(i) = ( dx * real(i-1,kind(dx)) )**3 - one
    end do
    !
  end subroutine SetExactSolution
  !
  subroutine AssembleA(this)
    !>
    class(gmres_solver_t) :: this
    !
    integer :: cnt, i
    !
  continue
    ! The 1st row
    cnt = 1
    i = 1
    this%ia(i) = cnt
    this%ja(cnt) = 1
    this%A(cnt) = one
    cnt = cnt + 1
    ! Only the interior nodes
    do i = 2, this%n-1
      this%ia(i) = cnt
      this%ja(cnt) = i - 1
      this%A(cnt) = this%idx2
      cnt = cnt + 1
      this%ja(cnt) = i
      this%A(cnt) = -two * this%idx2
      cnt = cnt + 1
      this%ja(cnt) = i + 1
      this%A(cnt) = this%idx2
      cnt = cnt + 1
    end do
    ! The last row
    i = this%n
    this%ia(i) = cnt
    this%ja(cnt) = this%n
    this%A(cnt) = one
    ! cnt = cnt + 1
    ! Store the real number of non-zero with 1 added
    this%ia(this%n+1) = cnt + 1
    !
  end subroutine AssembleA
!----------------------------------------------------------------------
! solve A x = b for x
! requires lapack routines dgelsy.
!----------------------------------------------------------------------
! m       gmres dimension
! n       dimension of x
! x       on input:  guess for x, can be 0
!         on exit:  solution x
! b       input b
! MatVec  performs y := A x, call MatVec(N,x, y)
! PCSolve  preconditioner, solve M x_out = x_in, call PCSolve(N,x)
! DotPrd  dot product, d = DotPrd(n,a,b)
! h       Hessian matrix,  size (m+1)*m
! v       Krylov subspace, size n*(m+1)
! res     on input: |Ax-b|/|b|<res;
!         on exit:  residual reached
! its     on input: max num its;
!         on exit:  number of its taken
! info    on input: if(info==1) print* residuals
!         on exit:  0 sucessful, 1 method breakdown, 2 max its
!----------------------------------------------------------------------
! Original by OpenPipeFlow
! Refer to the webpage: https://openpipeflow.org/index.php?title=File:GMRESm.f90
!----------------------------------------------------------------------

  subroutine Solve(slv)
    ! DGELSY computes the minimum-norm solution to a real linear least
    ! squares problem:
    ! minimize || A * X - B ||
    ! using a complete orthogonal factorization of A.  A is an M-by-N
    ! matrix which may be rank-deficient.
    ! For the example usage of DGELSY, check the webpage
    ! https://numericalalgorithmsgroup.github.io/LAPACK_Examples/examples/doc/dgelsy_example.html
    external :: DGELSY

    class(gmres_solver_t) :: slv
    !
    interface
      ! From SPARSKIT2
      subroutine AMUX(n, x, y, A, ja, ia)
        integer :: n
        integer, dimension(*) :: ja, ia
        real(8), dimension(*) :: A, x, y
      end subroutine AMUX
      !
    end interface
    !
    ! TODO: allocate the array w, z on the heap
    ! real(wp) :: tol,res_,stgn, w(n), z(n)
    ! real(wp) :: h_(m+1,m), y(m+1), p(m+1), work(4*m+1)
    ! integer :: imx, piv(m), rank, i
    real(wp) :: beta, beta_init, res_val, res_rel
    integer :: i, j, m, irst, rank
    logical :: print_info
    ! logical :: done
    real(wp), parameter :: rcond = 0.00001_wp
    !
  continue
    !
    m = slv%m
    if (slv%gmres_info==1) then
      print_info = .true.
    end if
    slv%its = 0
    ! With rst_max times restarting GMRES, it still fails to obtain an
    ! approx solution within the tolerance.
    ! Init to this value.
    slv%gmres_info = -1
    !
    ! loop_rst: do irst = 1, slv%rst_max
      ! Initialization
      ! write(*, *) slv%u
      call AMUX(slv%n, slv%u, slv%work_w, slv%A, slv%ja, slv%ia)
      ! write(*, *) slv%work_w
      slv%work_w(:) = slv%b - slv%work_w
      beta = norm2(slv%work_w)
      ! if (==1) then
        ! Save the initial residual
        beta_init = beta
      ! end if
      slv%kryl(:, 1) = slv%work_w / beta
      ! write(*, *) slv%kryl(:, 1)
      !
      slv%hess(:, :) = zero
      loop_m: do j = 1, m
        ! Right preconditioning
        ! z = M^{-1} v
        call slv%PCSolve(slv%kryl(:, j), slv%work_z)
        call AMUX(slv%n, slv%work_z, slv%work_w, slv%A, slv%ja, slv%ia)
        do i = 1, j
          slv%hess(i, j) = dot_product(slv%work_w, slv%kryl(:, i))
          slv%work_w(:) = slv%work_w - slv%hess(i, j) * slv%kryl(:, i)
        end do
        ! write(*, *) slv%work_w
        slv%hess(j+1, j) = norm2(slv%work_w)
        slv%kryl(:, j+1) = slv%work_w / slv%hess(j+1, j)
        ! write(*, *) slv%hess(:, j)
        ! write(*, *) slv%kryl(:, j+1)
        slv%work_y(1) = beta
        slv%work_y(2:j+1) = zero
        slv%hess_tmp(1:j+1, 1:j) = slv%hess(1:j+1, 1:j)
        ! write(*, *) slv%hess_tmp(1:j+1, 1:j)
        slv%dgelsy_piv = 0
        ! hess shall be altered by DGELSY. Thus we need an additional hess
        call DGELSY(j+1, j, 1, slv%hess_tmp(1:j+1, 1:j), j+1, slv%work_y(1:j+1), j+1, &
                    slv%dgelsy_piv, rcond, rank, slv%dgelsy_work, slv%dgelsy_lwork, i)
        if(i/=0) stop 'GMRESm: DGELSY'
        ! write(*, *) slv%work_y(1:j+1)
        ! Calculate the residual
        ! slv%work_y(:) = matmul(slv%hess(1:m+1, 1:m), slv%work_y(1:m))
        ! slv%work_y(1) = beta - slv%work_y(1)
        ! res_val = norm2(slv%work_y(1:m+1))
        res_val = abs(slv%work_y(j+1))
        res_rel = res_val / beta_init
        if(print_info) then
          write(*, 1) m, j, res_rel, res_val
        end if
        !
        if (res_val<=slv%atol .or. res_rel<=slv%rtol) then
          ! Obtain the approx solution within the tolerance
          slv%gmres_info = 0
          ! Exit
          exit loop_m
        else if (res_rel>slv%dtol) then
          ! Divergence
          slv%gmres_info = -2
          exit loop_m
        end if
        !
      end do loop_m
      !
      ! write(*, *) slv%work_y
      ! its = its + 1
      !
    ! end do loop_rst
    ! Update the solution
    if (slv%gmres_info==0 .or. slv%gmres_info==-2) then
      slv%work_w(:) = matmul(slv%kryl(:, 1:j), slv%work_y(1:j))
      call slv%PCSolve(slv%work_w, slv%work_z)
      slv%u(:) = slv%u + slv%work_z
      slv%its = j
    else if (slv%gmres_info==-1) then
      slv%work_w(:) = matmul(slv%kryl(:, 1:m), slv%work_y(1:m))
      call slv%PCSolve(slv%work_w, slv%work_z)
      slv%u(:) = slv%u + slv%work_z
      slv%its = m
    end if
    ! Impractical using 100000 iterations. I5 is enough for printing its
1 format('GMRES(m=', I0, '): restart no. ', I5, ', res(rel)=', F7.2, ', res(abs)=', ES14.6)
    !
  end subroutine Solve

  subroutine PCSolve(this, v_arr, w_arr)
    !
    class(gmres_solver_t) :: this
    !>
    real(wp), dimension(:), intent(in) :: v_arr
    !>
    real(wp), dimension(:), intent(out) :: w_arr
    !
  continue
    !
    if (this%gmres_pc_type==GMRES_No_PC) then
      ! No preconditioner
      w_arr(:) = v_arr
    else if (this%gmres_pc_type==GMRES_Jacobi) then
      ! Jacobi preconditioner
      w_arr(2:this%n-1) = this%idiag * v_arr(2:this%n-1)
      w_arr(1) = v_arr(1)
      w_arr(this%n) = v_arr(this%n)
    end if
    !
  end subroutine PCSolve

!   subroutine CalcResidual(this, r)
!     !
!     class(gmres_solver_t) :: this
!     !>
!     real(wp), dimension(:), intent(out) :: r
!     !
!     integer :: i
!     !
!   continue
!     ! r(0) and r(n) = 0, unchanged
!     do i = 2, this%n-1
!       r(i) = this%b(i) - ( this%u(i-1)+this%u(i+1)-two*this%u(i) ) * this%idx2
!     end do ! i
!     !
!   end subroutine CalcResidual

  subroutine ShowSaveResult(this)
    !
    class(gmres_solver_t) :: this
    !
    integer :: ierr, fn, i
    real(wp) :: norm
    !
  continue
    !
    open(newunit=fn, file=trim(this%res_fname), status="replace", &
         action="write")
    write(fn, *) "Vec Object: 1 MPI processes"
    write(fn, *) "  type: seq"
    do i = 1, this%n
      write(fn, 11) this%u(i)
    end do
    close(fn)
    !
    norm = norm2(this%u - this%u_exact)
    write(*, 100) norm, this%its
    !
11  format(ES23.15)
100 format('Norm of error ',ES12.4,' iterations ',I0)
    !
  end subroutine ShowSaveResult

!   subroutine DeconstructJacobiSolver(slv)
!     !
!     type(gmres_solver_t) :: slv
!     !
!     integer :: ierr
!     !
!   continue
!     !
!     !
!   end subroutine DeconstructJacobiSolver
!   !
  subroutine Free(this)
    !
    class(gmres_solver_t) :: this
    !
    integer :: ierr
    !
  continue
    !
    deallocate(this%u)
    deallocate(this%b)
    deallocate(this%u_exact)
    deallocate(this%A)
    deallocate(this%ja)
    deallocate(this%ia)
    deallocate(this%hess)
    deallocate(this%hess_tmp)
    deallocate(this%kryl)
    deallocate(this%work_w)
    deallocate(this%work_z)
    deallocate(this%work_y)
    !
  end subroutine Free
  !
end module gmres_solver_mod
