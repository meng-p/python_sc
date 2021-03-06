# from implicit import solver_im
from special_topic import solver_imsp
import matplotlib.pyplot as plt
# import sympy as sym
import numpy as np

dx = 0.01
L = 1
T = 0.1
dt = dx/2
mu = dt/(dx**2)


def u01(x):
    return np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)


# apply the implicit Euler scheme (sparce matrix version)
u, x, t = solver_imsp(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)

t = T
u_i = np.sin(np.pi*x)*np.exp(-(np.pi**2)*t) \
      + 0.1*np.sin(100*np.pi*x)*np.exp(-(100**2)*(np.pi**2)*t)

# compare the exact and numerical values
plt.plot(x, u, "-.", label="Finite difference approximation")
plt.plot(x, u_i, label="Exact solution")
plt.legend()
plt.show()
