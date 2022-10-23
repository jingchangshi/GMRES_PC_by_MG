import numpy as np

lsq_A = np.loadtxt('lsq.dat')
lsq_A = lsq_A.T
print(lsq_A)

lsq_b = np.zeros((5,))
#  lsq_b[0] = 8.9302855497458768
lsq_b[0] = 8.8652213975735581

lsq_x, lsq_res, A_rank, A_s = np.linalg.lstsq(lsq_A, lsq_b, rcond=None)
print(lsq_x)
print(lsq_res)
lsq_mul = np.matmul(lsq_A, lsq_x)
print(lsq_mul)
