import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('sjc')

#  ncell = 8
#  ncell = 16
ncell = 128
#  ncell = 1280
#  ncell = 10240
npts = ncell + 1
x_arr = np.linspace(0.0, 1.0, npts)
u_exact_arr = x_arr**3 - 1.0
#  gmres_ilu_arr = np.loadtxt('petsc_gmres_ilu.dat', skiprows=2)
petsc_gmres_none_arr = np.loadtxt('petsc_gmres_none.dat', skiprows=2)
petsc_gmres_jacobi_arr = np.loadtxt('petsc_gmres_jacobi.dat', skiprows=2)
#  petsc_gmres_sor_arr = np.loadtxt('petsc_gmres_sor.dat', skiprows=2)
petsc_gmres_mg_arr = np.loadtxt('petsc_gmres_mg.dat', skiprows=2)
#  petsc_gmres_shell_jacobi_arr = np.loadtxt('petsc_gmres_shell_jacobi.dat', skiprows=2)
#  jacobi_arr = np.loadtxt('jacobi.dat', skiprows=2)
#  multigrid_arr = np.loadtxt('multigrid.dat', skiprows=2)
#  gmres_no_pc_arr = np.loadtxt('gmres_no_pc.dat', skiprows=2)
#  gmres_jacobi_arr = np.loadtxt('gmres_jacobi.dat', skiprows=2)

fig=plt.figure()
ax=fig.gca()

marker_size = 2
plot_skip=2
#  plot_skip=16
ax.plot(x_arr, u_exact_arr, '-', label="Exact")
#  ax.plot(x_arr, gmres_ilu_arr, 's', ms=marker_size, label="GMRES-ILU")
ax.plot(x_arr[::plot_skip], petsc_gmres_none_arr[::plot_skip], '>', ms=marker_size, label="GMRES-None(PETSc)")
ax.plot(x_arr[::plot_skip], petsc_gmres_jacobi_arr[::plot_skip], '<', ms=marker_size, label="GMRES-Jacobi(PETSc)")
#  ax.plot(x_arr[::plot_skip], petsc_gmres_jacobi_arr[::plot_skip], 'v', ms=marker_size, label="GMRES-Sor(PETSc)")
ax.plot(x_arr[::plot_skip], petsc_gmres_mg_arr[::plot_skip], 's', ms=marker_size, label="GMRES-MG(PETSc)")
#  ax.plot(x_arr[::plot_skip], petsc_gmres_shell_jacobi_arr[::plot_skip], 's', ms=marker_size, label="GMRES-Jacobi(PETSc,Shell)")
#  ax.plot(x_arr[::plot_skip], jacobi_arr[::plot_skip], 'o', ms=marker_size, label="Jacobi")
#  ax.plot(x_arr[::plot_skip], multigrid_arr[::plot_skip], 'D', ms=marker_size, label="MG")
#  ax.plot(x_arr[::plot_skip], gmres_no_pc_arr[::plot_skip], 's', ms=marker_size, label="GMRES-NoPC")
#  ax.plot(x_arr[::plot_skip], gmres_jacobi_arr[::plot_skip], 'v', ms=marker_size, label="GMRES-Jacobi")
ax.legend(loc='lower right')

plt.savefig("result.png")
