module solver_mod
  !
  use dtype_mod
  !
  implicit none
  !
  private
  public :: solver_t
  !
  type, abstract :: solver_t
    !
    !
  contains
    !
    procedure(slv_GetN_t), deferred :: GetN
    procedure(slv_GetU_t), deferred :: GetU
    procedure(slv_SetU_t), deferred :: SetU
    procedure(slv_Solve_t), deferred :: Solve
    procedure(slv_SetB_t), deferred :: SetB
    procedure(slv_CalcResidual_t), deferred :: CalcResidual
    procedure(slv_ShowSaveResult_t), deferred :: ShowSaveResult
    !
  end type solver_t
  !
  interface
    !
    function slv_GetN_t(this) result(n)
      import solver_t
      class(solver_t) :: this
      integer :: n
    end function slv_GetN_t
    !
    subroutine slv_GetU_t(this, u)
      import solver_t, wp
      class(solver_t) :: this
      real(wp), dimension(:), intent(out) :: u
    end subroutine slv_GetU_t
    !
    subroutine slv_SetU_t(this, u)
      import solver_t, wp
      class(solver_t) :: this
      real(wp), dimension(:), intent(in) :: u
    end subroutine slv_SetU_t
    !
    subroutine slv_Solve_t(this)
      import solver_t
      class(solver_t) :: this
    end subroutine slv_Solve_t
    !
    subroutine slv_SetB_t(this, r)
      import solver_t, wp
      class(solver_t) :: this
      real(wp), dimension(:), intent(in), optional :: r
    end subroutine slv_SetB_t
    !
    subroutine slv_CalcResidual_t(this, r)
      import solver_t, wp
      class(solver_t) :: this
      real(wp), dimension(:), intent(out) :: r
    end subroutine slv_CalcResidual_t
    !
    subroutine slv_ShowSaveResult_t(this)
      import solver_t
      class(solver_t) :: this
    end subroutine slv_ShowSaveResult_t
    !
  end interface

end module solver_mod
