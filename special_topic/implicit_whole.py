import numpy as np


def solver_im2(alpha, In, L, T, dt, dx, mu):
    Nt = int(round(T/float(dt)))
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)
    t = np.linspace(0, T, Nt+1)
    u = np.zeros((Nt+1, Nx+1))
    A = np.zeros((Nx+1, Nx+1))
    b = np.zeros(Nx+1)
    A[0, 0] = A[Nx, Nx] = 1
    b[0] = b[Nx] = 0
    for i in range(1, Nx):
        A[i, i-1] = -mu
        A[i, i+1] = -mu
        A[i, i] = 1 + 2*mu
    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u[0, i] = In(x[i])
    import scipy.linalg
    for n in range(1, Nt+1):
        for i in range(1, Nx):
            b[i] = u[n-1, i]
        u[n, :] = scipy.linalg.solve(A, b)

    return u, x, t
