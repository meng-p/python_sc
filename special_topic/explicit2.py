import numpy as np

# forward in time, central in space


def solver_em2(alpha, In, L, T, dt, dx, mu):
    # import time
    # t0 = time.perf_counter()
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, T, Nt+1)   # Mesh points in time
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)
    u = np.zeros((Nt+1, Nx+1))  # store the results at every time level
    # u[0,] = np.zeros(Nx+1)

    for i in range(0, Nx+1):
        u[0, i] = In(x[i])
    for m in range(1, Nt+1):  # apply the explicit Euler scheme
        for j in range(1, Nx):
            u[m, j] = u[m-1, j] \
             + mu*(u[m-1, j+1] - 2*u[m-1, j] + u[m-1, j-1])

    # t1 = time.perf_counter()
    return u, x, t
