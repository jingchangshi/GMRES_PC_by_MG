module jacobi_solver_mod
  !
  use dtype_mod
  !
  implicit none
  !
  private
  public :: jacobi_solver_t
  !
  type :: jacobi_solver_t
    !
    integer :: n
    real(wp), dimension(:), allocatable :: u
    real(wp), dimension(:), allocatable :: b
    real(wp), dimension(:), allocatable :: u_nxt
    real(wp), dimension(:), allocatable :: u_exact
    integer :: maxits
    !
    real(wp) :: idx2, idiag
    integer :: its
    character(len=64) :: res_fname
    !
  contains
    !
    procedure :: SetExactSolution
    procedure :: GetU, SetU, Solve, SetB, CalcResidual, ShowSaveResult, GetN, Free
    procedure :: SolveAsPreconditioner
    !
  end type jacobi_solver_t
  !
  interface jacobi_solver_t
    module procedure ConstructJacobiSolver
  end interface jacobi_solver_t
  !
contains
  !
  function ConstructJacobiSolver(n, maxits) result(slv)
    !>
    integer, intent(in) :: n
    !>
    integer, intent(in) :: maxits
    !>
    type(jacobi_solver_t) :: slv
    !
    integer :: i
    real(wp) :: dx
    !
  continue
    !
    slv%n = n
    slv%maxits = maxits
    !
    dx = one / real(slv%n-1, kind(dx))
    !
    allocate(slv%u(1:slv%n), source=zero)
    allocate(slv%b(1:slv%n), source=zero)
    allocate(slv%u_nxt(1:slv%n), source=zero)
    allocate(slv%u_exact(1:slv%n), source=zero)
    !
    call slv%SetB()
    call slv%SetExactSolution()
    !
    slv%idx2 = one / dx**2
    slv%idiag = -one / ( two * slv%idx2 )
    !
    slv%res_fname = "jacobi.dat"
    !
  end function ConstructJacobiSolver
  !
  subroutine SetB(this, r)
    !>
    class(jacobi_solver_t) :: this
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
    class(jacobi_solver_t) :: this
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
  subroutine GetU(this, u)
    !>
    class(jacobi_solver_t) :: this
    !
    real(wp), dimension(:), intent(out) :: u
    !
  continue
    !
    u(:) = this%u
    !
  end subroutine GetU
  !
  subroutine SetU(this, u)
    !>
    class(jacobi_solver_t) :: this
    !
    real(wp), dimension(:), intent(in) :: u
    !
  continue
    !
    this%u(:) = u
    !
  end subroutine SetU
  !
  subroutine Solve(this)
    !
    class(jacobi_solver_t) :: this
    !
    integer :: i, j
    !
  continue
    ! the first node
    this%u(1) = this%b(1)
    ! the last node
    this%u(this%n) = this%b(this%n)
    do i = 1, this%maxits
      do j = 2, this%n-1
        this%u_nxt(j) = ( this%b(j) - (this%u(j-1)+this%u(j+1))*this%idx2 ) * &
                        this%idiag
      end do
      this%u(2:this%n-1) = this%u_nxt(2:this%n-1)
      !
    end do ! i
    !
    this%its = i - 1
    !
  end subroutine Solve

  subroutine SolveAsPreconditioner(this, r, du)
    !
    class(jacobi_solver_t) :: this
    !
    real(wp), dimension(:), intent(in) :: r
    !
    real(wp), dimension(:), intent(out) :: du
    !
    integer :: n
    !
  continue
    !
    n = size(du,dim=1)
    du(2:n-1) = this%idiag * r(2:n-1)
    du(1) = r(1)
    du(n) = r(n)
    ! this%its = i - 1
    !
  end subroutine SolveAsPreconditioner

  subroutine CalcResidual(this, r)
    !
    class(jacobi_solver_t) :: this
    !>
    real(wp), dimension(:), intent(out) :: r
    !
    integer :: i
    !
  continue
    ! r(0) and r(n) = 0, unchanged
    do i = 2, this%n-1
      r(i) = this%b(i) - ( this%u(i-1)+this%u(i+1)-two*this%u(i) ) * this%idx2
    end do ! i
    !
  end subroutine CalcResidual

  subroutine ShowSaveResult(this)
    !
    class(jacobi_solver_t) :: this
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

  function GetN(this) result(n)
    !
    class(jacobi_solver_t) :: this
    !
    integer :: n
    !
  continue
    !
    n = this%n
    !
  end function GetN

  subroutine DeconstructJacobiSolver(slv)
    !
    type(jacobi_solver_t) :: slv
    !
    integer :: ierr
    !
  continue
    !
    !
  end subroutine DeconstructJacobiSolver
  !
  subroutine Free(this)
    !
    class(jacobi_solver_t) :: this
    !
    integer :: ierr
    !
  continue
    !
    deallocate(this%u)
    deallocate(this%b)
    deallocate(this%u_nxt)
    deallocate(this%u_exact)
    !
  end subroutine Free
  !
end module jacobi_solver_mod
