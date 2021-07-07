import numpy as np

n = 8
A = np.zeros((n+1, n+1))
u = np.zeros((n+1,))
b = np.zeros((n+1,))

dx = (1.0-0.0)/n

for i in range(1,n):
    A[i, i-1] = 1 / dx**2
    A[i, i  ] =-2 / dx**2
    A[i, i+1] = 1 / dx**2
i = 0
A[i, 0] = 1
i = n
A[i, n] = 1

for i in range(1,n):
    b[i] = i*dx * 6.0
b[0] =-1.0
b[n] = 0.0

u = np.linalg.solve(A, b)
