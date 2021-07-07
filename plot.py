import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('sjc')

ncell = 10240
npts = ncell + 1
x_arr = np.linspace(0.0, 1.0, npts)
u_exact_arr = x_arr**3 - 1.0
#  gmres_ilu_arr = np.loadtxt('gmres_ilu.dat', skiprows=2)
#  gmres_jacobi_arr = np.loadtxt('gmres_jacobi.dat', skiprows=2)
#  gmres_matrdef_jacobi_arr = np.loadtxt('gmres_matrdef_jacobi.dat', skiprows=2)
#  gmres_userdef_jacobi_arr = np.loadtxt('gmres_userdef_jacobi.dat', skiprows=2)
jacobi_arr = np.loadtxt('jacobi.dat', skiprows=2)
multigrid_arr = np.loadtxt('multigrid.dat', skiprows=2)

fig=plt.figure()
ax=fig.gca()

marker_size = 2

ax.plot(x_arr, u_exact_arr, '-', label="Exact")
#  ax.plot(x_arr, gmres_ilu_arr, 's', ms=marker_size, label="GMRES-ILU")
#  ax.plot(x_arr, gmres_jacobi_arr, 's', ms=marker_size, label="GMRES-Jacobi")
#  ax.plot(x_arr, gmres_matrdef_jacobi_arr, 'v', ms=marker_size, label="GMRES-Jacobi(MatrDef)")
#  ax.plot(x_arr, gmres_userdef_jacobi_arr, 's', ms=marker_size, label="GMRES-Jacobi(UserDef)")
ax.plot(x_arr, jacobi_arr, 'o', ms=marker_size, label="Jacobi")
ax.plot(x_arr, multigrid_arr, 'v', ms=marker_size, label="MG")
ax.legend()

plt.savefig("result.png")
