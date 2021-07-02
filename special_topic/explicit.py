import numpy as np

# forward in time, central in space


def solver_em(alpha, In, L, T, dt, dx, mu):
    import time
    t0 = time.perf_counter()  # measure the CPU time
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, T, Nt+1)   # grid points in time
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)  # grid points in space
    u = np.zeros(Nx+1)  # solution at new time level (unknown)
    u_n = np.zeros(Nx+1)  # u at previous time level

    for j in range(0, Nx+1):
        u_n[j] = In(x[j])  # set the initial condition u(x,0)=In(x)

    # use an explicit loop for solving the finite difference scheme
    for n in range(0, Nt):  # compute the solution at the internal points
        for j in range(1, Nx):
            u[j] = u_n[j] + mu*(u_n[j+1] - 2*u_n[j] + u_n[j-1])  # apply EES
        u[0] = 0
        u[Nx] = 0  # boundary conditions
        u_n, u = u, u_n  # update from the last time level

    t1 = time.perf_counter()
    return u_n, x, t, t1-t0  # u_n holds latest u
