program main
  !
! #include <petsc/finclude/petscksp.h>
!   use petscksp
  use dtype_mod
  use solver_mod, only: solver_t
  use petsc_solver_mod, only: petsc_solver_t
  use jacobi_solver_mod, only: jacobi_solver_t
  use multigrid_solver_mod, only: multigrid_solver_t
  !
  implicit none
  !
  class(solver_t), allocatable :: solver
  !>
  integer :: n
  !>
  real(wp) :: rtol
  !>
  real(wp) :: atol
  !>
  real(wp) :: dtol
  !>
  integer :: maxits
  !>
  integer :: gmres_pc_type
  !>
  integer :: solver_type
  !>
  character(len=64) :: level_maxits_input_file
  !>
  integer :: argc
  !>
  character(len=64), dimension(:), allocatable :: argv
  !
  integer :: i, j, fn
  integer :: nlevel, ncycle, mg_slv_type
  integer, dimension(:), allocatable :: level_maxits
  character(len=64) :: str
  !
!  Note: Any user-defined Fortran routines MUST be declared as external.
  ! external MatrDef_Jacobi_ShellPCSetUp, MatrDef_Jacobi_ShellPCApply, &
  !          MatrDef_Jacobi_ShellPCDestroy
  ! external UserDef_Jacobi_ShellPCSetUp, UserDef_Jacobi_ShellPCApply, &
  !          UserDef_Jacobi_ShellPCDestroy

  continue
  ! Default values
  ! Number of cells
  n = 8
  rtol = 1E-7_wp
  atol = 1E-14_wp
  dtol = 1E2_wp
  maxits = 1000
  gmres_pc_type = PC_Jacobi
  solver_type = PETSc
  ! Update based on the cmd line argument
  argc = command_argument_count()
  if (mod(argc,2)/=0) then
    write(*, *) "Invalid number of arguments!"
    stop
  end if
  !
  allocate(argv(1:argc))
  do i = 1, argc
    call get_command_argument(i, str)
    argv(i) = str
    ! write(*, *) argv
  end do
  !
  do i = 0, argc/2-1
    j = i * 2 + 1
    select case (argv(j))
    case ("-n")
      ! Number of cells
      read(argv(j+1), *) n
    case ("-rtol")
      read(argv(j+1), *) rtol
    case ("-atol")
      read(argv(j+1), *) atol
    case ("-dtol")
      read(argv(j+1), *) dtol
    case ("-maxits")
      read(argv(j+1), *) maxits
    case ("-gmres_pc_type")
      read(argv(j+1), *) gmres_pc_type
    case ("-solver_type")
      read(argv(j+1), *) solver_type
    case ("-nlevel")
      read(argv(j+1), *) nlevel
    case ("-ncycle")
      read(argv(j+1), *) ncycle
    case ("-level_maxits_input_file")
      read(argv(j+1), *) level_maxits_input_file
    case default
      write(*, 101) argv(j)
      stop
    end select
  end do
101 format("Invalid argument: ", A)
  !
  deallocate(argv)
  ! n is the number of points [0, 1, ..., n]
  n = n + 1
  !
  if (solver_type==Multigrid) then
    mg_slv_type = Jacobi
    ! allocate(level_solvers(1:nlevel), source=mg_slv_type)
    allocate(level_maxits(1:nlevel), source=0)
    open(newunit=fn, file=trim(level_maxits_input_file), status="old", action="read")
    do j = 1, nlevel
      read(fn, *) level_maxits(j)
    end do
    close(fn)
  end if ! 
  !
  ! Show the simulation info
  !
  write(*, 11) "n", n
  write(*, 11) "solver_type", solver_type
  if (solver_type==Multigrid) then
    write(*, 11) "nlevel", nlevel
    write(*, 11) "ncycle", ncycle
    write(*, '(A)', advance='no') "level_maxits is ["
    do j = 1, nlevel-1
      write(*, 21, advance='no') level_maxits(j)
    end do
    write(*, 22) level_maxits(nlevel)
  else if (solver_type==Jacobi) then
    write(*, 11) "maxits", maxits
  end if
11 format(A, " is ", I0)
21 format(I0, ", ")
22 format(I0, "]")
  !
  if (solver_type==PETSc) then
    allocate(solver, source = &
             petsc_solver_t(n, rtol, atol, dtol, maxits, gmres_pc_type))
  else if (solver_type==Jacobi) then
    allocate(solver, source = &
             jacobi_solver_t(n, maxits))
  else if (solver_type==Multigrid) then
    allocate(solver, source = &
             multigrid_solver_t( &
               nlevel, mg_slv_type, level_maxits, n, ncycle))
  else
    write(*, *) "solver_type invalid!"
    stop
  end if
  !
  call solver%Solve()
  call solver%ShowSaveResult()
  !
  deallocate(solver)
  !
  ! dx = (one - zero) / real(n, kind(dx))
  ! ! In the serial mode, Istart = 0, Iend = n
  ! !
  ! solver_type = PETSc
  ! call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
  !                         '-solver_type',solver_type,flg,ierr)
  ! gmres_pc_type = PC_Jacobi
  ! call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
  !                         '-gmres_pc_type',gmres_pc_type,flg,ierr)
  ! rtol = 1E-7_wp
  ! call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
  !                         '-rtol',rtol,flg,ierr)
  ! atol = 1E-14_wp
  ! call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
  !                         '-atol',atol,flg,ierr)
  ! dtol = 1E+2_wp
  ! call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
  !                         '-dtol',dtol,flg,ierr)
  ! maxits = 100
  ! call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
  !                         '-maxits',maxits,flg,ierr)
  ! !
  ! col_val(-1) = one / dx**2
  ! col_val( 0) =-two / dx**2
  ! col_val( 1) = one / dx**2
  ! !
  ! call AssembleVector(u_exact, n, ierr)
  ! ! call VecView(b, PETSC_VIEWER_STDOUT_SELF, ierr)
  ! ! call VecView(u, PETSC_VIEWER_STDOUT_SELF, ierr)
  ! if (solver_type==PETSc) then
  !   !
  !   ! PETSc solver
  !   !

  !   if (any(gmres_pc_type == PC_UserDefs) .or. &
  !       any(gmres_pc_type == PC_MatrDefs)) then
  !     call PCSetType(pc, PCSHELL, ierr)
  !     if (gmres_pc_type==PC_MatrDef_Jacobi) then
  !       write(*, *) "Using PC_MatrDef_Jacobi"
  !       res_fname = "gmres_matrdef_jacobi.dat"
  !       call PCShellSetApply(pc, MatrDef_Jacobi_ShellPCApply, ierr)
  !       call PCShellSetSetUp(pc, MatrDef_Jacobi_ShellPCSetUp, ierr)
  !       call PCShellSetDestroy(pc, MatrDef_Jacobi_ShellPCDestroy, ierr)
  !     else if (gmres_pc_type==PC_UserDef_Jacobi) then
  !       write(*, *) "Using PC_UserDef_Jacobi"
  !       res_fname = "gmres_userdef_jacobi.dat"
  !       call PCShellSetApply(pc, UserDef_Jacobi_ShellPCApply, ierr)
  !       call PCShellSetSetUp(pc, UserDef_Jacobi_ShellPCSetUp, ierr)
  !       call PCShellSetDestroy(pc, UserDef_Jacobi_ShellPCDestroy, ierr)
  !     else
  !       write(*,*) "gmres_pc_type invalid!"
  !       stop
  !     end if
  !   else
  !   endif

  !   ! call VecView(u, PETSC_VIEWER_STDOUT_SELF, ierr)
  !   ! call VecView(u_exact, PETSC_VIEWER_STDOUT_SELF, ierr)

  !   !
  ! else if (solver_type==Jacobi) then
  !   !
  !   ! Jacobi method
  !   !
  !   write(*, *) "Using Jacobi method to solve Ax=b"
  !   res_fname = "jacobi.dat"
  !   ! Set up the matrix and vector for the Jacobi method
  !   call VecDuplicate(u, tmp_u, ierr)
  !   !
  !   call VecDuplicate(u, diag, ierr)
  !   call MatGetDiagonal(A, diag, ierr)
  !   call VecReciprocal(diag, ierr)
  !   !
  !   call MatDuplicate(A, MAT_COPY_VALUES, LU, ierr)
  !   do II = Istart, Iend-1
  !     call MatSetValue(LU, II, II, zero, INSERT_VALUES, ierr)
  !   end do
  !   call MatAssemblyBegin(LU, MAT_FINAL_ASSEMBLY, ierr)
  !   call MatAssemblyEnd(LU, MAT_FINAL_ASSEMBLY, ierr)
  !   ! call MatView(LU, PETSC_VIEWER_STDOUT_SELF, ierr)
  !   !
  !   !
  !   its = i-1
  !   !
  ! else
  !   write(*, *) "solver_type invalid!"
  !   stop
  ! end if

end program
!
! The Jacobi preconditioner with an explicit form
!
! subroutine MatrDef_Jacobi_ShellPCSetUp(pc, ierr)
!   use petsc_data_mod
!   PetscErrorCode :: ierr
!   PC :: pc
!   Mat :: pmat
! continue
!   call PCGetOperators(pc, PETSC_NULL_MAT, pmat, ierr)
!   call MatCreateVecs(pmat, diag, PETSC_NULL_VEC, ierr)
!   call MatGetDiagonal(pmat, diag, ierr)
!   call VecReciprocal(diag, ierr)
! end subroutine MatrDef_Jacobi_ShellPCSetUp

! subroutine MatrDef_Jacobi_ShellPCApply(pc, x, y, ierr)
!   use petsc_data_mod
!   PetscErrorCode :: ierr
!   PC :: pc
!   Vec :: x,y
! continue
!   call VecPointwiseMult(y, x, diag, ierr)
! end subroutine MatrDef_Jacobi_ShellPCApply

! subroutine MatrDef_Jacobi_ShellPCDestroy(pc, ierr)
!   use petsc_data_mod
!   PetscErrorCode :: ierr
!   PC :: pc
! continue
!   call VecDestroy(diag, ierr)
! end subroutine MatrDef_Jacobi_ShellPCDestroy
! !
! ! User defined Jacobi preconditioner without an explicit matrix form
! !
! subroutine UserDef_Jacobi_ShellPCSetUp(pc, ierr)
!   use petsc_data_mod
!   PetscErrorCode :: ierr
!   PC :: pc
!   Mat :: pmat
!   PetscInt :: II, n
!   integer, parameter :: wp = 8
!   real(wp), parameter :: zero = 0.0_wp, one = 1.0_wp, six = 6.0_wp
!   PetscReal :: dx
! continue
!   ! Set up the matrix and vector for the Jacobi method
!   call PCGetOperators(pc, PETSC_NULL_MAT, pmat, ierr)
!   call MatCreateVecs(pmat, diag, PETSC_NULL_VEC, ierr)
!   call MatGetDiagonal(pmat, diag, ierr)
!   call VecReciprocal(diag, ierr)
!   !
!   call VecGetSize(diag, n, ierr)
!   dx = (one - zero) / real(n-1, kind(dx))
!   call VecDuplicate(diag, b, ierr)
!   do II = 1, n-2
!     call VecSetValue(b, II, six*dx*real(II,kind(dx)), &
!                      INSERT_VALUES, ierr)
!   end do
!   call VecSetValue(b, 0, -one, INSERT_VALUES, ierr)
!   call VecSetValue(b, n-1, zero, INSERT_VALUES, ierr)
!   !
!   call MatDuplicate(pmat, MAT_COPY_VALUES, LU, ierr)
!   do II = 0, n-1
!     call MatSetValue(LU, II, II, zero, INSERT_VALUES, ierr)
!   end do
!   call MatAssemblyBegin(LU, MAT_FINAL_ASSEMBLY, ierr)
!   call MatAssemblyEnd(LU, MAT_FINAL_ASSEMBLY, ierr)
!   call MatView(LU, PETSC_VIEWER_STDOUT_SELF, ierr)
!   call VecView(diag, PETSC_VIEWER_STDOUT_SELF, ierr)
!   call VecView(b, PETSC_VIEWER_STDOUT_SELF, ierr)
! end subroutine UserDef_Jacobi_ShellPCSetUp

! subroutine UserDef_Jacobi_ShellPCApply(pc, x, y, ierr)
!   use petsc_data_mod
!   PetscErrorCode :: ierr
!   PC :: pc
!   Vec :: x,y
!   integer :: i
!   integer, parameter :: wp = 8
!   integer, parameter :: maxits = 1
!   real(wp), parameter :: one = 1.0_wp
! continue

!   do i = 1, maxits
!     call MatMult(LU, x, y, ierr)
!     call VecAYPX(y, -one, b, ierr)
!     call VecPointwiseMult(y, diag, y, ierr)
!   end do ! i

! end subroutine UserDef_Jacobi_ShellPCApply

! subroutine UserDef_Jacobi_ShellPCDestroy(pc, ierr)
!   use petsc_data_mod
!   PetscErrorCode :: ierr
!   PC :: pc
! continue
!   call VecDestroy(diag, ierr)
!   call VecDestroy(b, ierr)
!   call MatDestroy(LU, ierr)
! end subroutine UserDef_Jacobi_ShellPCDestroy


