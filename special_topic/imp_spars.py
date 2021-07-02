import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg.dsolve.linsolve import spsolve


def solver_imsp(alpha, In, L, T, dt, dx, mu):
    Nt = int(round(T/float(dt)))
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)  # mesh points in space
    t = np.linspace(0, T, Nt+1)  # mesh points in time
    u = np.zeros(Nx+1)
    u_n = np.zeros(Nx+1)

    main = np.zeros(Nx+1)  # the sparse matrix
    lower = np.zeros(Nx)
    upper = np.zeros(Nx)
    b = np.zeros(Nx+1)  # right-hand side vector
    # the matrix for implicit Euler scheme
    main[:] = 1 + 2*mu  # main diagonal
    lower[:] = -mu
    upper[:] = -mu
    main[0] = 1  # insert boundary conditions
    main[Nx] = 1
    lower[0] = 0
    lower[-1] = 0
    upper[0] = 0
    upper[-1] = 0
    A = diags(diagonals=[main, lower, upper], offsets=[0, -1, 1],
              shape=(Nx+1, Nx+1), format='csr')
    for j in range(0, Nx+1):
        u_n[j] = In(x[j])
    for m in range(0, Nt):
        for j in range(1, Nx):
            b[j] = u_n[j]
        b[0] = b[Nx] = 0
        # solve the sparce linear system using the function 'spsolve'
        u[:] = spsolve(A, b)
        u_n[:] = u
    return u_n, x, t
