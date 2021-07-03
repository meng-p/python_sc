from special_topic import solver_cn1
import matplotlib.pyplot as plt
# import sympy as sym
import numpy as np

L = 1
T = 0.1
# u01 = lambda x: np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)


def u01(x):
    return np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)


# mu = 0.5
E = np.zeros(3)
dxx = [0.01, 0.05, 0.1]

for i in range(0, 3):
    dx = dxx[i]  # Nx = int(round(L/dx))
    # dt = mu*dx**2
    dt = dx/2
    mu = dt/(dx**2)
    u, x, t = solver_cn1(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
    t = T
    # x = L - dx
    u_e = np.sin(np.pi*x)*np.exp(-(np.pi**2)*t) \
        + 0.1*np.sin(100*np.pi*x)*np.exp(-(100**2)*(np.pi**2)*t)

    # E[i] = abs(u[-2] - u_e)
    E[i] = np.linalg.norm(u-u_e, np.inf)

# dtt = [0.5*0.01**2, 0.5*0.05**2, 0.5*0.1**2]
# check the convergence rate in time using loglog plot
dtt = [0.01/2, 0.05/2, 0.1/2]
dxxx = [0.01**2, 0.05**2, 0.1**2]
plt.loglog(dtt, E, label="convergence rate")
plt.loglog(dtt, np.square(dtt), '--', label="slope=2")
plt.loglog(dtt, dtt, '--', label="slope=1")
plt.xlabel(r"$\Delta$ t")
plt.ylabel("Error")
plt.legend()
plt.show()
