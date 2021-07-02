import numpy as np


def solver_cn1(alpha, In, L, T, dt, dx, mu):
    Nt = int(round(T/float(dt)))
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)  # grid points in space
    t = np.linspace(0, T, Nt+1)  # grid points in time
    u = np.zeros(Nx+1)  # unknown u at the new time level
    u_n = np.zeros(Nx+1)  # u at the previous time level
    A = np.zeros((Nx+1, Nx+1))  # matrix to store the results
    b = np.zeros(Nx+1)  # the right-hand side of the linear system
    for i in range(1, Nx):
        A[i, i-1] = -0.5*mu
        A[i, i+1] = -0.5*mu
        A[i, i] = 1 + mu
    A[0, 0] = A[Nx, Nx] = 1
    # Set initial condition u(x,0) = I(x)
    for i in range(0, Nx+1):
        u_n[i] = In(x[i])
    import scipy.linalg
    for n in range(0, Nt):  # compute the vector b
        for i in range(1, Nx):
            b[i] = u_n[i] + mu*0.5*(u_n[i+1]-2*u_n[i]+u_n[i-1])
        b[0] = b[Nx] = 0
        # solve the linear system using function 'solve'
        u[:] = scipy.linalg.solve(A, b)
        u_n[:] = u  # update u_n before the next step
    return u_n, x, t
