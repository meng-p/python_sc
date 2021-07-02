from special_topic import solver_em, solver_emv
import numpy as np


def test_solver_ex():
    dx = 0.01
    L = 1
    T = 0.1
    Nx = int(round(L/dx))
    a = 1
    dt = 0.5*(dx**2)
    mu = a*dt/(dx**2)  # mu=0.5
    x = np.linspace(0, L, Nx+1)

    def u01(x):
        return np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)

    u, x, t, cpu = solver_em(In=u01, alpha=a, L=L, T=T, dx=dx, dt=dt, mu=mu)
    t = T
    u_e = np.sin(np.pi*x)*np.exp(-(np.pi**2)*t) \
        + 0.1*np.sin(100*np.pi*x)*np.exp(-(100**2)*(np.pi**2)*t)
    diff = abs(u_e - u).max()
    tol = 1E-4
    assert diff < tol, 'max diff solver_FE_simple: %g' % diff

    u, x, t, cpu = solver_emv(In=u01, alpha=a, L=L, T=T, dx=dx, dt=dt, mu=mu)
    diff = abs(u_e - u).max()
    tol = 1E-4
    assert diff < tol, 'max diff solver_FE, scalar: %g' % diff
