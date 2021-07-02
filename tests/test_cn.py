from special_topic import solver_cn1
import numpy as np
from special_topic.cn_sparse import solver_cnsp


def test_solver_cn():  # verification
    dx = 0.01
    L = 1
    T = 0.1
    Nx = int(round(L/dx))
    a = 1
    dt = dx/2
    mu = a*dt/(dx**2)  # mu=0.5
    x = np.linspace(0, L, Nx+1)

    def u01(x):
        return np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)

    u, x, t = solver_cn1(In=u01, alpha=a, L=L, T=T, dx=dx, dt=dt, mu=mu)
    t = T
    # define the exact solution
    u_e = np.sin(np.pi*x)*np.exp(-(np.pi**2)*t) \
        + 0.1*np.sin(100*np.pi*x)*np.exp(-(100**2)*(np.pi**2)*t)
    diff1 = abs(u_e - u).max()
    tol = 1E-4
    assert diff1 < tol, 'max diff: %g' % diff1

    u, x, t = solver_cnsp(In=u01, alpha=a, L=L, T=T, dx=dx, dt=dt, mu=mu)
    diff = abs(u_e - u).max()
    tol = 1E-4
    assert diff < tol, 'max diff: %g' % diff

    # check that the vectoerised version has the same results
    assert diff1 - diff < tol
