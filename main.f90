program main
  !
  use dtype_mod
  use petsc_solver_mod, only: petsc_solver_t
  use jacobi_solver_mod, only: jacobi_solver_t
  use multigrid_solver_mod, only: multigrid_solver_t
  use gmres_solver_mod, only: gmres_solver_t
  !
  implicit none
  !
  type(jacobi_solver_t), allocatable :: jacobi_solver
  type(multigrid_solver_t), allocatable :: multigrid_solver
  type(petsc_solver_t), allocatable :: petsc_solver
  type(gmres_solver_t), allocatable :: gmres_solver
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
  !> GMRES dimension
  integer :: m
  !>
  integer :: argc
  !>
  character(len=64), dimension(:), allocatable :: argv
  !
  integer :: i, j, fn
  integer :: nlevel, ncycle, mg_slv_type, gmres_info
  integer, dimension(:), allocatable :: level_maxits
  character(len=64) :: str
  !
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
  m = 8
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
    case ("-m")
      read(argv(j+1), *) m
    case default
      ! write(*, 101) argv(j)
      ! stop
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
  if (solver_type==Jacobi) then
    allocate(jacobi_solver, source=jacobi_solver_t(n, maxits))
  else if (solver_type==Multigrid) then
    allocate(multigrid_solver, source = &
             multigrid_solver_t(nlevel, mg_slv_type, level_maxits, n, ncycle))
  else if (solver_type==PETSc) then
    allocate(petsc_solver, source = &
             petsc_solver_t(n, rtol, atol, dtol, maxits, gmres_pc_type))
  else if (solver_type==GMRES) then
    gmres_info = 1
    allocate(gmres_solver, source = &
             gmres_solver_t(m, n, gmres_pc_type, maxits, rtol, atol, dtol, gmres_info))
  else
    write(*, *) "solver_type invalid!"
    stop
  end if
  !
  if (solver_type==Jacobi) then
    call jacobi_solver%Solve()
    call jacobi_solver%ShowSaveResult()
    deallocate(jacobi_solver)
  else if (solver_type==PETSc) then
    call petsc_solver%Solve()
    call petsc_solver%ShowSaveResult()
    deallocate(petsc_solver)
  else if (solver_type==Multigrid) then
    call multigrid_solver%Solve()
    call multigrid_solver%ShowSaveResult()
    deallocate(multigrid_solver)
  else if (solver_type==GMRES) then
    call gmres_solver%Solve()
    call gmres_solver%ShowSaveResult()
    deallocate(gmres_solver)
  end if
  !
end program
