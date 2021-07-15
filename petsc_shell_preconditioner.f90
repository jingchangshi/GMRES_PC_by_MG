module petsc_shell_data_mod
#include <petsc/finclude/petscksp.h>
  use petscksp
  ! use dtype_mod

  use jacobi_solver_mod, only: jacobi_solver_t
  ! use multigrid_solver_mod, only: multigrid_solver_t

  type(jacobi_solver_t), allocatable :: jacobi_solver
  ! type(multigrid_solver_t), allocatable :: multigrid_solver

  ! real(wp), dimension(:), allocatable :: pc_shell_sol

end module petsc_shell_data_mod
!
! Jacobi solver as the preconditioner
!
subroutine PC_Shell_Jacobi_SetUp(pc, ierr)
  !
  use petsc_shell_data_mod
  !
  PC :: pc
  PetscErrorCode :: ierr
  !
  Mat :: pmat
  integer :: n
  integer, parameter :: maxits = 1
  ! type(pc_shell_ctx_t), pointer :: ctx_ptr
  ! PetscFortranAddr :: ctx
  !
continue
  !
  call PCGetOperators(pc, PETSC_NULL_MAT, pmat, ierr)
  call MatGetSize(pmat, n, PETSC_NULL_INTEGER, ierr)
  allocate(jacobi_solver, source=jacobi_solver_t(n, maxits))
  !
  ! allocate(pc_shell_sol(1:n), source=zero)
  ! call PCShellGetContext(pc, ctx, ierr)
  ! write(*, *) ctx
  !
end subroutine PC_Shell_Jacobi_SetUp

subroutine PC_Shell_Jacobi_Apply(pc, x, y, ierr)
  !
  use petsc_shell_data_mod
  !
  PC :: pc
  Vec :: x
  Vec :: y
  PetscErrorCode :: ierr
  !
  integer, parameter :: wp = 8
  real(wp), dimension(:), pointer :: x_ptr, y_ptr
  !
continue
  ! Set up the matrix and vector for the Jacobi method
  call VecGetArrayReadF90(x, x_ptr, ierr)
  call VecGetArrayF90(y, y_ptr, ierr)
  call jacobi_solver%SolveAsPreconditioner(x_ptr, y_ptr)
  call VecRestoreArrayReadF90(x, x_ptr, ierr)
  call VecRestoreArrayF90(y, y_ptr, ierr)
  !
end subroutine PC_Shell_Jacobi_Apply

subroutine PC_Shell_Jacobi_Destroy(pc, ierr)
  !
  use petsc_shell_data_mod
  !
  PC :: pc
  PetscErrorCode :: ierr
  !
continue
  !
  deallocate(jacobi_solver)
  !
end subroutine PC_Shell_Jacobi_Destroy
!
! Multigrid solver as the preconditioner
!
! subroutine PC_Shell_MG_SetUp(pc, ierr)
!   !
!   use petsc_shell_data_mod
!   !
!   PC :: pc
!   PetscErrorCode :: ierr
!   !
!   Mat :: pmat
!   integer :: n, nlevel, ncycle, mg_slv_type
!   integer, dimension(:), allocatable :: level_maxits
!   integer :: fn, j
!   character(len=64) :: level_maxits_input_file
!   !
! continue
!   !
!   nlevel = 3
!   !
!   mg_slv_type = Jacobi
!   !
!   level_maxits_input_file = "level_maxits.in"
!   allocate(level_maxits(1:nlevel), source=0)
!   open(newunit=fn, file=trim(level_maxits_input_file), status="old", action="read")
!   do j = 1, nlevel
!     read(fn, *) level_maxits(j)
!   end do
!   close(fn)
!   !
!   call PCGetOperators(pc, PETSC_NULL_MAT, pmat, ierr)
!   call MatGetSize(pmat, n, PETSC_NULL_INTEGER, ierr)
!   !
!   ncycle = 1
!   !
!   allocate(multigrid_solver, source = &
!            multigrid_solver_t(nlevel, mg_slv_type, level_maxits, n, ncycle))
!   !
!   deallocate(level_maxits)
!   !
! end subroutine PC_Shell_MG_SetUp

! subroutine PC_Shell_MG_Apply(pc, x, y, ierr)
!   !
!   use petsc_shell_data_mod
!   !
!   PC :: pc
!   Vec :: x
!   Vec :: y
!   PetscErrorCode :: ierr
!   !
!   integer, parameter :: wp = 8
!   real(wp), dimension(:), pointer :: x_ptr, y_ptr
!   !
! continue
!   ! Set up the matrix and vector for the Jacobi method
!   call VecGetArrayReadF90(x, x_ptr, ierr)
!   call VecGetArrayF90(y, y_ptr, ierr)
!   call multigrid_solver%SolveAsPreconditioner(x_ptr, y_ptr)
!   call VecRestoreArrayReadF90(x, x_ptr, ierr)
!   call VecRestoreArrayF90(y, y_ptr, ierr)
!   !
! end subroutine PC_Shell_MG_Apply

! subroutine PC_Shell_MG_Destroy(pc, ierr)
!   !
!   use petsc_shell_data_mod
!   !
!   PC :: pc
!   PetscErrorCode :: ierr
!   !
! continue
!   !
!   deallocate(multigrid_solver)
!   !
! end subroutine PC_Shell_MG_Destroy

