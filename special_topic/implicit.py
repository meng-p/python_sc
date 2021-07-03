import numpy as np
from scipy.linalg import solve


def solver_im(In, L, T, dt, dx, mu):
    Nt = int(round(T/float(dt)))
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)
    t = np.linspace(0, T, Nt+1)
    u = np.zeros(Nx+1)
    u_n = np.zeros(Nx+1)
    A = np.zeros((Nx+1, Nx+1))  # construct the linear system
    b = np.zeros(Nx+1)  # the right-hand side of the linear system
    for j in range(1, Nx):
        A[j, j-1] = -mu
        A[j, j+1] = -mu
        A[j, j] = 1 + 2*mu
    A[0, 0] = A[Nx, Nx] = 1
    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u_n[i] = In(x[i])
    for m in range(0, Nt):
        for j in range(1, Nx):
            b[j] = u_n[j]
        b[0] = b[Nx] = 0
        u[:] = solve(A, b)
        u_n[:] = u  # switch variables before the next step
    return u_n, x, t
