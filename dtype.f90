module dtype_mod

  public

  integer, parameter :: wp = 8
  real(wp), parameter :: zero = 0.0_wp, one = 1.0_wp, two = 2.0_wp, &
                         half = 0.5_wp, one4 = 0.25_wp, six = 6.0_wp, &
                         eps12 = 1E-12_wp

  integer, parameter :: PETSc = 1, Jacobi = 2, Multigrid = 3
  integer, parameter :: PC_Jacobi = 1, PC_SOR = 2, PC_ILU = 3, PC_MG = 4, &
                        PC_Shell_Jacobi = 11
  ! integer, parameter :: PC_Shell_MG = 12
  ! integer, parameter :: PC_UserDefs(1:1) = [PC_UserDef_Jacobi]
  !
  type :: mat1r_t
    real(wp), dimension(:), allocatable :: v
  end type mat1r_t

end module dtype_mod
