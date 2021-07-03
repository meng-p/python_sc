from special_topic import solver_cn1
import matplotlib.pyplot as plt
# import sympy as sym
import numpy as np

dx = 0.01  # step size in space
L = 1
T = 0.1
dt = dx/2
mu = dt/(dx**2)  # mu=0.5


def u01(x):  # initial condition
    return np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)


# apply the Crank-Nicolson scheme
u, x, t = solver_cn1(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)

t = T
u_e = np.sin(np.pi*x)*np.exp(-(np.pi**2)*t) \
      + 0.1*np.sin(100*np.pi*x)*np.exp(-(100**2)*(np.pi**2)*t)


# compare the exact solution and the numerical solution
plt.plot(x, u, "-.", label="Finite difference approximation")
plt.plot(x, u_e, label="Exact solution")
plt.legend()
plt.show()
