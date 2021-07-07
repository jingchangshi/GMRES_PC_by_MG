module multigrid_solver_mod
  !
  use dtype_mod
  use solver_mod, only: solver_t
  use jacobi_solver_mod, only: jacobi_solver_t
  !
  implicit none
  !
  private
  public :: multigrid_solver_t
  !
  type :: mg_slv_t
    class(solver_t), allocatable :: ptr
  end type mg_slv_t
  !
  type, extends(solver_t) :: multigrid_solver_t
    !
    type(mg_slv_t), dimension(:), allocatable :: solvers
    integer :: nlevel
    ! integer :: n
    integer :: maxcycle
    !
    type(mat1r_t), dimension(:), allocatable :: r
    type(mat1r_t), dimension(:), allocatable :: u
    type(mat1r_t), dimension(:), allocatable :: z
    !
    character(len=64) :: res_fname
    !
  contains
    !
    procedure :: Prolongate, Restrict
    procedure :: GetU, SetU, Solve, SetB, CalcResidual, ShowSaveResult, GetN, Free
    !
  end type multigrid_solver_t
  !
  interface multigrid_solver_t
    module procedure ConstructMultigridSolver
  end interface multigrid_solver_t
  !
contains
  !
  function ConstructMultigridSolver( &
      nlevel, slv_type, level_maxits, &
      nnode_finest, maxcycle) result(slv)
    !
    use jacobi_solver_mod, only: jacobi_solver_t
    !>
    integer, intent(in) :: nlevel
    !>
    integer, intent(in) :: slv_type
    !>
    integer, dimension(:), intent(in) :: level_maxits
    !>
    integer, intent(in) :: nnode_finest
    !>
    integer, intent(in) :: maxcycle
    !>
    type(multigrid_solver_t) :: slv
    !
    integer :: i, n
    real(wp) :: dx
    !
  continue
    !
    slv%nlevel = nlevel
    slv%maxcycle = maxcycle
    select case (slv_type)
    case (Jacobi)
      allocate(slv%solvers(1:slv%nlevel))
    case (PETSc)
    case default
      write(*, *) "Invalid solver type!"
      stop
    end select
    !
    allocate(slv%r(1:slv%nlevel))
    allocate(slv%u(1:slv%nlevel))
    allocate(slv%z(1:slv%nlevel))
    !
    do i = 1, slv%nlevel
      !
      n = (nnode_finest-1) / 2**(i-1) + 1
      !
      select case (slv_type)
      case (Jacobi)
        slv%solvers(i)%ptr = jacobi_solver_t(n, level_maxits(i))
      case (PETSc)
      case default
        write(*, *) "Invalid solver type!"
        stop
      end select
      !
      allocate(slv%r(i)%v(1:n), source=zero)
      allocate(slv%u(i)%v(1:n), source=zero)
      allocate(slv%z(i)%v(1:n), source=zero)
      !
    end do ! i
    !
    slv%res_fname = "multigrid.dat"
    !
  end function ConstructMultigridSolver
  !
  subroutine Prolongate(this, l, h, arr_l, arr_h)
    !
    class(multigrid_solver_t) :: this
    !>
    integer, intent(in) :: l
    !>
    integer, intent(in) :: h
    !>
    real(wp), dimension(:), intent(in) :: arr_l
    !>
    real(wp), dimension(:), intent(out) :: arr_h
    !
    integer :: i, j, nl, nh
    !
  continue
    !
    nl = this%solvers(l)%ptr%GetN()
    nh = this%solvers(h)%ptr%GetN()
    !
    do i = 2, nl-1
      j = i*2
      arr_h(j-1) = arr_l(i)
      arr_h(j) = (arr_l(i) + arr_l(i+1)) * half
    end do ! i
    arr_h(2) = (arr_l(1) + arr_l(2)) * half
    arr_h(nh-1) = (arr_l(nl-1) + arr_l(nl)) * half
    !
  end subroutine Prolongate
  !
  subroutine Restrict(this, h, l, arr_h, arr_l)
    !
    class(multigrid_solver_t) :: this
    !>
    integer, intent(in) :: h
    !>
    integer, intent(in) :: l
    !>
    real(wp), dimension(:), intent(in) :: arr_h
    !>
    real(wp), dimension(:), intent(out) :: arr_l
    !
    integer :: i, j, nl, nh
    !
  continue
    !
    nl = this%solvers(l)%ptr%GetN()
    nh = this%solvers(h)%ptr%GetN()
    !
    do i = 2, nl-1
      j = (i-1)*2+1
      arr_l(i) = one4 * (arr_h(j-1) + two*arr_h(j) + arr_h(j+1))
    end do ! i
    !
  end subroutine Restrict

  subroutine GetU(this, u)
    !
    class(multigrid_solver_t) :: this
    !>
    real(wp), dimension(:), intent(out) :: u
    !
  continue
    !
  end subroutine GetU

  subroutine SetU(this, u)
    !
    class(multigrid_solver_t) :: this
    !>
    real(wp), dimension(:), intent(in) :: u
    !
  continue
    !
  end subroutine SetU

  subroutine Solve(this)
    !
    class(multigrid_solver_t) :: this
    !>
    ! real(wp), dimension(:), intent(in), optional :: r
    !
    integer :: ic, il, it, i, j
    !
  continue
    !
    do ic = 1, this%maxcycle
      ! Loop all levels down to the coarest level
      do il = 1, this%nlevel-1
        ! Pre-smooth
        call this%solvers(il)%ptr%Solve()
        ! Restrict the residual
        call this%solvers(il)%ptr%CalcResidual(this%r(il)%v)
        call this%Restrict(il, il+1, this%r(il)%v, this%r(il+1)%v)
        ! this%solvers(il+1)%ptr%b(:) = this%r(il+1)%v
        call this%solvers(il+1)%ptr%SetB(this%r(il+1)%v)
        call this%solvers(il+1)%ptr%SetU(this%z(il+1)%v)
      end do ! il
      ! Solve on the coarest level
      call this%solvers(il)%ptr%Solve()
      ! Loop all levels up to the finest level
      do il = this%nlevel-1, 1, -1
        ! Prolongate
        ! call this%solvers(il+1)%ptr%CalcResidual()
        call this%solvers(il+1)%ptr%GetU(this%u(il+1)%v)
        ! To save memory, use r for e
        call this%Prolongate(il+1, il, this%u(il+1)%v, this%r(il)%v)
        call this%solvers(il)%ptr%GetU(this%u(il)%v)
        this%u(il)%v = this%u(il)%v + this%r(il)%v
        call this%solvers(il)%ptr%SetU(this%u(il)%v)
        ! Post-smooth
        call this%solvers(il)%ptr%Solve()
      end do ! il
      !
    end do ! ic
    !
  end subroutine Solve

  subroutine SetB(this, r)
    !
    class(multigrid_solver_t) :: this
    !>
    real(wp), dimension(:), intent(in), optional :: r
    !
  continue
    !
  end subroutine SetB

  subroutine CalcResidual(this, r)
    !
    class(multigrid_solver_t) :: this
    !>
    real(wp), dimension(:), intent(out) :: r
    !
  continue
    !
  end subroutine CalcResidual

  subroutine ShowSaveResult(this)
    !
    class(multigrid_solver_t) :: this
    !
    integer :: ierr, fn, i, n
    real(wp) :: norm
    !
  continue
    !
    open(newunit=fn, file=trim(this%res_fname), status="replace", &
         action="write")
    write(fn, *) "Vec Object: 1 MPI processes"
    write(fn, *) "  type: seq"
    !
    call this%solvers(1)%ptr%GetU(this%u(1)%v)
    n = this%solvers(1)%ptr%GetN()
    do i = 1, n
      write(fn, 11) this%u(1)%v(i)
    end do
    close(fn)
    !
    ! norm = norm2(this%u - this%u_exact)
    ! write(*, 100) norm, this%its
    !
11  format(ES23.15)
100 format('Norm of error ',ES12.4,' iterations ',I0)
    !
  end subroutine ShowSaveResult

  function GetN(this) result(n)
    !
    class(multigrid_solver_t) :: this
    !
    integer :: n
    !
  continue
    !
    ! n = this%n
    !
  end function GetN

  subroutine DeconstructMultigridSolver(slv)
    !
    type(multigrid_solver_t) :: slv
    !
    integer :: ierr
    !
  continue
    !
    !
  end subroutine DeconstructMultigridSolver
  !
  subroutine Free(this)
    !
    class(multigrid_solver_t) :: this
    !
    integer :: ierr
    !
  continue
    !
    !
  end subroutine Free
  !
end module multigrid_solver_mod
