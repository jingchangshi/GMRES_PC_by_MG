FC=mpiifort
FFLAGS=-fpp -g -O0
PETSc_DIR=/home/jcshi/NFS_Project/NFS/contrib_intel/petsc_dbg
# FC=mpif90
# FFLAGS=-cpp -g -O0
# PETSc_DIR=/home/jcshi/NFS_Project/NFS/contrib_gnu/petsc_dbg
PETSc_FLAGS=-I${PETSc_DIR}/include -L${PETSc_DIR}/lib -lpetsc -lflapack -lfblas -lHYPRE -ldl
solver: main.f90 dtype.o petsc_solver.o jacobi_solver.o multigrid_solver.o
	${FC} ${FFLAGS} -o solver main.f90 dtype.o petsc_solver.o jacobi_solver.o multigrid_solver.o ${PETSc_FLAGS}
dtype.o: dtype.f90
	${FC} ${FFLAGS} -c dtype.f90
# solver.o: solver.f90 dtype.o
#   ${FC} ${FFLAGS} -c solver.f90 ${PETSc_FLAGS}
jacobi_solver.o: jacobi_solver.f90
	${FC} ${FFLAGS} -c jacobi_solver.f90 ${PETSc_FLAGS}
petsc_solver.o: petsc_solver.f90
	${FC} ${FFLAGS} -c petsc_solver.f90 ${PETSc_FLAGS}
multigrid_solver.o: multigrid_solver.f90 jacobi_solver.o petsc_solver.o
	${FC} ${FFLAGS} -c multigrid_solver.f90
clean:
	rm solver *.o *.mod

run01:
	./solver -n 10240 -solver_type 1 -gmres_pc_type 3 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 1000
run02:
	./solver -n 100 -solver_type 1 -gmres_pc_type 1 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 1000
run03:
	./solver -n 100 -solver_type 1 -gmres_pc_type 2 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 1000
run11:
	./solver -n 100 -solver_type 1 -gmres_pc_type 11 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 1000
run21:
	./solver -n 100 -solver_type 1 -gmres_pc_type 21 -rtol 1E-7 -atol 1E-14 -dtol 1E2 -maxits 10000
run31:
	./solver -n 10240 -solver_type 2 -maxits 20000
run41:
	./solver -n 10240 -solver_type 3 -level_maxits_input_file level_maxits.in -nlevel 6 -ncycle 1000

