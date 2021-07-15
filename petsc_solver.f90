module petsc_solver_mod
  !
#include <petsc/finclude/petscksp.h>
  use petscksp
  !
  use dtype_mod
  !
  implicit none
  !
  private
  public :: petsc_solver_t
  !
  PetscInt, parameter :: nnz_per_row = 3
  !
  type :: petsc_solver_t
    !
    PetscInt :: n
    Mat :: A
    Vec :: u, b
    KSP :: ksp
    PC :: pc
    Vec :: u_exact
    !
    PetscReal :: idx2
    !
    PetscReal :: rtol, atol, dtol
    PetscInt :: maxits
    PetscInt :: gmres_pc_type
    !
    PetscInt :: its
    PetscViewer :: viewer
    character(len=64) :: res_fname
    !
    contains
    !
    procedure :: AssembleMatrix, AssembleVector, SetExactSolution
    procedure :: GetU, SetU, Solve, SetB, CalcResidual, ShowSaveResult, GetN, Free
    final :: DeconstructPetscSolver
    !
  end type petsc_solver_t
  !
  interface petsc_solver_t
    module procedure ConstructPetscSolver
  end interface petsc_solver_t
  !
!  Note: Any user-defined Fortran routines MUST be declared as external.
  external PC_Shell_Jacobi_Apply, PC_Shell_Jacobi_SetUp, PC_Shell_Jacobi_Destroy
  ! external PC_Shell_MG_Apply, PC_Shell_MG_SetUp, PC_Shell_MG_Destroy
  !
contains
  !
  function ConstructPetscSolver(n, rtol, atol, dtol, maxits, gmres_pc_type &
            ) result(slv)
    !
    ! use petsc_shell_data_mod, only: petsc_pc_shell_ctx
    !>
    integer, intent(in) :: n
    !>
    real(wp), intent(in) :: rtol
    !>
    real(wp), intent(in) :: atol
    !>
    real(wp), intent(in) :: dtol
    !>
    integer, intent(in) :: maxits
    !>
    integer, intent(in) :: gmres_pc_type
    !>
    type(petsc_solver_t) :: slv
    !
    PetscErrorCode :: ierr
    !
  continue
    !
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    if (ierr /= 0) then
      write(*, *) 'Unable to initialize PETSc'
      stop
    endif
    !
    slv%n = n
    slv%rtol = rtol
    slv%atol = atol
    slv%dtol = dtol
    slv%maxits = maxits
    slv%gmres_pc_type = gmres_pc_type
    !
    write(*, 101) slv%n
101 format("Using PETSc solver with the number of grid points: ", I0)
    !
    call slv%AssembleMatrix()
    ! call MatView(A, PETSC_VIEWER_STDOUT_SELF, ierr)
    ! Assemble the solution vector u
    call slv%AssembleVector(slv%u)
    ! Set the initial solution value
    call VecSet(slv%u, zero, ierr)
    call slv%AssembleVector(slv%b)
    call slv%SetB()
    !
    call slv%AssembleVector(slv%u_exact)
    call slv%SetExactSolution()
    !
    call KSPCreate(PETSC_COMM_WORLD, slv%ksp, ierr)
    call KSPSetType(slv%ksp, KSPGMRES, ierr)
    call KSPSetOperators(slv%ksp, slv%A, slv%A, ierr)
    call KSPSetTolerances(slv%ksp, slv%rtol, slv%atol, slv%dtol, slv%maxits, ierr)
    !
    call KSPGetPC(slv%ksp, slv%pc, ierr)
    if (slv%gmres_pc_type==PC_Jacobi) then
      write(*, *) "Using PCJACOBI"
      slv%res_fname = "gmres_jacobi.dat"
      call PCSetType(slv%pc, PCJACOBI, ierr)
    else if (slv%gmres_pc_type==PC_SOR) then
      write(*, *) "Using PCSOR"
      slv%res_fname = "gmres_sor.dat"
      call PCSetType(slv%pc, PCSOR, ierr)
    else if (slv%gmres_pc_type==PC_ILU) then
      write(*, *) "Using PCILU"
      slv%res_fname = "gmres_ilu.dat"
      call PCSetType(slv%pc, PCILU, ierr)
    else if (slv%gmres_pc_type==PC_MG) then
      write(*, *) "Using PCMG"
      slv%res_fname = "gmres_mg.dat"
      call PCSetType(slv%pc, PCMG, ierr)
      ! call PCMGSetType(slv%pc, PC_MG_MULTIPLICATIVE, ierr)
      ! call PCMGSetCycleType(slv%pc, PC_MG_CYCLE_V, ierr)
    else if (slv%gmres_pc_type==PC_Shell_Jacobi) then
      write(*, *) "Using PC_Shell_Jacobi"
      slv%res_fname = "gmres_shell_jacobi.dat"
      call PCSetType(slv%pc, PCSHELL, ierr)
      ! slv%pc_shell_ctx = slv%b
      ! call VecDuplicate(slv%b,slv%pc_shell_ctx,ierr)
      ! call VecCopy(slv%b,slv%pc_shell_ctx,ierr)
      ! allocate(slv%pc_shell_ctx_ptr)
      ! petsc_pc_shell_ctx%n = 107
      ! slv%pc_shell_ctx_ptr => petsc_pc_shell_ctx
      ! call PCShellSetContext(slv%pc, slv%pc_shell_ctx_ptr, ierr)
      ! call PCShellSetContext(slv%pc, slv%pc_shell_ctx, ierr)
      call PCShellSetSetUp(slv%pc, PC_Shell_Jacobi_SetUp, ierr)
      call PCShellSetApply(slv%pc, PC_Shell_Jacobi_Apply, ierr)
      call PCShellSetDestroy(slv%pc, PC_Shell_Jacobi_Destroy, ierr)
    else
      write(*, *) "gmres_pc_type invalid!"
      stop
    end if
    call KSPSetFromOptions(slv%ksp, ierr)
    !
    ! call MatView(slv%A, PETSC_VIEWER_STDOUT_SELF, ierr)
    ! call VecView(slv%u, PETSC_VIEWER_STDOUT_SELF, ierr)
    ! call VecView(slv%b, PETSC_VIEWER_STDOUT_SELF, ierr)
    !
  end function ConstructPetscSolver
  !
  subroutine AssembleMatrix(this)
    !>
    class(petsc_solver_t) :: this
    !
    PetscErrorCode :: ierr
    !
    PetscInt :: II, col_idx(-1:1)
    PetscReal :: dx, col_val(-1:1)
    !
  continue
    !
    call MatCreate(PETSC_COMM_WORLD, this%A, ierr)
    call MatSetSizes(this%A, PETSC_DECIDE, PETSC_DECIDE, this%n, this%n, ierr)
    call MatSetType(this%A, MATAIJ, ierr)
    call MatSetFromOptions(this%A, ierr)
    ! For now only serial running.
    ! call MatMPIAIJSetPreallocation(A,3,PETSC_NULL_INTEGER,five,PETSC_NULL_INTEGER,ierr)
    call MatSeqAIJSetPreallocation(this%A, nnz_per_row, PETSC_NULL_INTEGER, ierr)
    !
    dx = one / real(this%n-1, kind(dx))
    this%idx2 = one / dx**2
    col_val(-1) = this%idx2
    col_val( 0) =-two * this%idx2
    col_val( 1) = this%idx2
    ! call MatGetOwnershipRange(A, Istart, Iend, ierr)
    ! Only the interior nodes
    do II = 1, this%n-2
      col_idx(-1) = II-1
      col_idx( 0) = II
      col_idx( 1) = II+1
      call MatSetValues(this%A, 1, II, 3, col_idx, col_val, INSERT_VALUES, ierr)
    end do
    ! The 1st row
    call MatSetValues(this%A, 1, 0, 1, 0, one, INSERT_VALUES, ierr)
    ! The last row
    call MatSetValues(this%A, 1, this%n-1, 1, this%n-1, one, INSERT_VALUES, ierr)
    !
    call MatAssemblyBegin(this%A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(this%A, MAT_FINAL_ASSEMBLY, ierr)
    !
  end subroutine AssembleMatrix

  subroutine AssembleVector(this, x)
    !>
    class(petsc_solver_t) :: this
    !>
    Vec :: x
    !
    PetscErrorCode :: ierr
    !
  continue
    !
    call VecCreate(PETSC_COMM_WORLD, x, ierr)
    call VecSetSizes(x, PETSC_DECIDE, this%n, ierr)
    call VecSetFromOptions(x, ierr)
    !
  end subroutine AssembleVector

  subroutine SetB(this, r)
    !>
    class(petsc_solver_t) :: this
    !
    real(wp), dimension(:), intent(in), optional :: r
    !
    PetscErrorCode :: ierr
    PetscInt :: II
    PetscReal :: dx
    !
  continue
    !
    if (present(r)) then
      do II = 0, this%n-1
        call VecSetValue(this%b, II, r(II+1), INSERT_VALUES, ierr)
      end do
    else
      dx = one / real(this%n-1, kind(dx))
      do II = 1, this%n-2
        call VecSetValue(this%b, II, six*dx*real(II,kind(dx)), INSERT_VALUES, ierr)
      end do
      call VecSetValue(this%b, 0, -one, INSERT_VALUES, ierr)
      call VecSetValue(this%b, this%n-1, zero, INSERT_VALUES, ierr)
    end if
    !
    call VecAssemblyBegin(this%b, ierr)
    call VecAssemblyEnd(this%b, ierr)
    !
  end subroutine SetB

  subroutine SetExactSolution(this)
    !>
    class(petsc_solver_t) :: this
    !
    PetscErrorCode :: ierr
    PetscInt :: II
    PetscReal :: dx
    !
  continue
    dx = one / real(this%n-1, kind(dx))
    ! Analytic solution
    do II = 0, this%n-1
      call VecSetValue(this%u_exact, II, (dx*real(II,kind(dx)))**3-one, &
                       INSERT_VALUES, ierr)
    end do
    call VecAssemblyBegin(this%u_exact, ierr)
    call VecAssemblyEnd(this%u_exact, ierr)
    !
  end subroutine SetExactSolution

  subroutine GetU(this, u)
    !>
    class(petsc_solver_t) :: this
    !
    real(wp), dimension(:), intent(out) :: u
    !
    PetscErrorCode :: ierr
    PetscInt :: II
    PetscScalar, dimension(:), pointer :: u_ptr
    !
  continue
    !
    call VecGetArrayReadF90(this%u, u_ptr, ierr)
    do II = 0, this%n-1
      u(II+1) = u_ptr(II)
    end do
    call VecRestoreArrayReadF90(this%u, u_ptr, ierr)
    !
  end subroutine GetU

  subroutine SetU(this, u)
    !>
    class(petsc_solver_t) :: this
    !
    real(wp), dimension(:), intent(in) :: u
    !
    PetscErrorCode :: ierr
    PetscInt :: II
    PetscScalar, dimension(:), pointer :: u_ptr
    !
  continue
    !
    call VecGetArrayF90(this%u, u_ptr, ierr)
    do II = 0, this%n-1
      u_ptr(II) = u(II+1)
    end do
    call VecRestoreArrayF90(this%u, u_ptr, ierr)
    !
  end subroutine SetU

  subroutine Solve(this)
    !
    class(petsc_solver_t) :: this
    !
    integer :: i
    !
    PetscErrorCode :: ierr
    !
  continue
    !
    call KSPSolve(this%ksp, this%b, this%u, ierr)
    call KSPGetIterationNumber(this%ksp, this%its, ierr)
    !
  end subroutine Solve

  subroutine CalcResidual(this, r)
    !
    class(petsc_solver_t) :: this
    !>
    real(wp), dimension(:), intent(out) :: r
    !
    PetscErrorCode :: ierr
    PetscScalar, dimension(:), pointer :: b_ptr, u_ptr
    integer :: i
    !
  continue
    !
    call VecGetArrayReadF90(this%u, u_ptr, ierr)
    call VecGetArrayReadF90(this%b, b_ptr, ierr)
    ! r(0) and r(n) = 0, unchanged
    ! The index starts from 0 but r starts from 1
    do i = 1, this%n-2
      r(i+1) = b_ptr(i) - ( u_ptr(i-1)+u_ptr(i+1)-two*u_ptr(i) ) * this%idx2
    end do ! i
    !
    call VecRestoreArrayReadF90(this%u, u_ptr, ierr)
    call VecRestoreArrayReadF90(this%b, b_ptr, ierr)
    !
  end subroutine CalcResidual

  subroutine ShowSaveResult(this)
    !
    class(petsc_solver_t) :: this
    !
    PetscErrorCode :: ierr
    PetscReal :: norm
    !
  continue
    !
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD, trim(this%res_fname), &
                              this%viewer, ierr)
    call VecView(this%u, this%viewer, ierr)
    !
    call VecAXPY(this%u, -one, this%u_exact, ierr)
    call VecNorm(this%u, NORM_2, norm, ierr)
    write(*, 100) norm, this%its
    !
100 format('Norm of error ',ES12.4,' iterations ',I0)
    !
  end subroutine ShowSaveResult

  function GetN(this) result(n)
    !
    class(petsc_solver_t) :: this
    !
    integer :: n
    !
  continue
    !
    n = this%n
    !
  end function GetN

  subroutine DeconstructPetscSolver(slv)
    !
    type(petsc_solver_t) :: slv
    !
    PetscErrorCode :: ierr
    !
  continue
    !
    ! call KSPDestroy(slv%ksp, ierr)
    ! call VecDestroy(slv%u, ierr)
    ! call VecDestroy(slv%b, ierr)
    ! call VecDestroy(slv%u_exact, ierr)
    ! call MatDestroy(slv%A, ierr)

    ! call PetscFinalize(ierr)
    !
  end subroutine DeconstructPetscSolver
  !
  subroutine Free(this)
    !
    class(petsc_solver_t) :: this
    !
    PetscErrorCode :: ierr
    !
  continue
    !
    call KSPDestroy(this%ksp, ierr)
    call VecDestroy(this%u, ierr)
    call VecDestroy(this%b, ierr)
    call VecDestroy(this%u_exact, ierr)
    call MatDestroy(this%A, ierr)

    call PetscFinalize(ierr)
    !
  end subroutine Free
  !
end module petsc_solver_mod
