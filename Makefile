# INTEL
FC=mpiifort
FFLAGS=-fpp -g -O0
PETSc_DIR=/home/jcshi/NFS_Project/NFS/contrib_intel/petsc_dbg
PETSc_FLAGS=-I${PETSc_DIR}/include -L${PETSc_DIR}/lib -lpetsc -lflapack -lfblas -lHYPRE -ldl
SPARSKIT2_DIR=/home/jcshi/NFS_Project/NFS/contrib_intel/SPARSKIT2_dbg
SPARSKIT2_FLAGS=-L${SPARSKIT2_DIR}/lib -lskit
LAPACK_DIR=/home/jcshi/NFS_Project/NFS/contrib_intel/LAPACK
LAPACK_FLAGS=-L${LAPACK_DIR}/lib64 -llapack -lblas
# GNU
# FC=mpif90
# FFLAGS=-cpp -g -O0
# PETSc_DIR=/home/jcshi/NFS_Project/NFS/contrib_gnu/petsc_dbg
# PETSc_FLAGS=-I${PETSc_DIR}/include -L${PETSc_DIR}/lib -lpetsc -lflapack -lfblas -lHYPRE -ldl
# SPARSKIT2_DIR=/home/jcshi/NFS_Project/NFS/contrib_gnu/SPARSKIT2_dbg
# SPARSKIT2_FLAGS=-L${SPARSKIT2_DIR}/lib -lskit
# LAPACK_DIR=/home/jcshi/NFS_Project/NFS/contrib_gnu/LAPACK
# LAPACK_FLAGS=-L${LAPACK_DIR}/lib64 -llapack -lblas
solver: main.f90 dtype.o petsc_solver.o jacobi_solver.o multigrid_solver.o petsc_shell_preconditioner.o gmres_solver.o
	${FC} ${FFLAGS} -o solver main.f90 dtype.o petsc_solver.o jacobi_solver.o multigrid_solver.o petsc_shell_preconditioner.o gmres_solver.o ${PETSc_FLAGS} ${SPARSKIT2_FLAGS} ${LAPACK_FLAGS}
dtype.o: dtype.f90
	${FC} ${FFLAGS} -c dtype.f90
# solver.o: solver.f90 dtype.o
#   ${FC} ${FFLAGS} -c solver.f90 ${PETSc_FLAGS}
jacobi_solver.o: jacobi_solver.f90
	${FC} ${FFLAGS} -c jacobi_solver.f90 ${PETSc_FLAGS}
petsc_solver.o: petsc_solver.f90 petsc_shell_preconditioner.o
	${FC} ${FFLAGS} -c petsc_solver.f90 ${PETSc_FLAGS}
multigrid_solver.o: multigrid_solver.f90 jacobi_solver.o
	${FC} ${FFLAGS} -c multigrid_solver.f90
petsc_shell_preconditioner.o: petsc_shell_preconditioner.f90 jacobi_solver.o multigrid_solver.o
	${FC} ${FFLAGS} -c petsc_shell_preconditioner.f90 ${PETSc_FLAGS}
gmres_solver.o: gmres_solver.f90
	${FC} ${FFLAGS} -c gmres_solver.f90 ${SPARSKIT2_FLAGS} ${LAPACK_FLAGS}
clean:
	rm solver *.o *.mod

run11:
	./solver -n 1280 -solver_type 1 -gmres_pc_type 1 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 1000
run12:
	./solver -n 1280 -solver_type 1 -gmres_pc_type 2 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 1000
run13:
	./solver -n 1280 -solver_type 1 -gmres_pc_type 3 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 1000
run14:
	./solver -n 1280 -solver_type 1 -gmres_pc_type 4 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 1000
run111:
	./solver -n 1280 -solver_type 1 -gmres_pc_type 11 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 1000
# run112:
#   ./solver -n 1280 -solver_type 1 -gmres_pc_type 12 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 10000
run21:
	./solver -n 1280 -solver_type 2 -maxits 100
run31:
	./solver -n 1280 -solver_type 3 -level_maxits_input_file level_maxits.in -nlevel 6 -ncycle 500
run41:
	./solver -n 1280 -solver_type 4 -gmres_pc_type 1 -m 100 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 100

