from special_topic import solver_em
import matplotlib.pyplot as plt
# import sympy as sym
import numpy as np

L = 1
T = 0.1


def u01(x):
    return np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)
# u01 = lambda x: np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)


# mu = 0.5
E = np.zeros(3)
dxx = [0.01, 0.05, 0.1]

for i in range(0, 3):
    dx = dxx[i]  # Nx = int(round(L/dx))
    # dt = mu*dx**2
    dt = 0.5*(dx**2)
    # mu = a*dt/(dx**2)
    mu = dt/(dx**2)
    u, x, t, cpu = solver_em(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
    t = T
    # x = L - dx
    u_e = np.sin(np.pi*x)*np.exp(-(np.pi**2)*t) \
        + 0.1*np.sin(100*np.pi*x)*np.exp(-(100**2)*(np.pi**2)*t)

    # E[i] = abs(u[-2] - u_e)
    E[i] = np.linalg.norm(u-u_e, np.inf)

dxxx = [0.01**2, 0.05**2, 0.1**2]
# the convergence rate of the explicit method in space
plt.loglog(dxx, E, label="convergence rate")
plt.loglog(dxx, dxxx, '--', label="slope=2")
plt.loglog(dxx, dxx, '--', label="slope=1")
plt.xlabel(r"$\Delta$ x")
plt.ylabel("Error")
plt.legend()
plt.show()
