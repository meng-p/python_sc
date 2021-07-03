from special_topic import solver_emv, solver_cnsp, solver_imsp
# import sympy as sym
import matplotlib.pyplot as plt
import numpy as np

dx = 0.1  # set the step size in space
# dx = 0.01
L = 1
# L = np.pi
T = 0.1  # a given final time
# T = 3
dt = 0.5*(dx**2)  # the step size in time depends on dx
mu = 0.5  # mu=0.5


def u01(x):  # initial condition
    return np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)
    # return np.sin(x)
# u01 = lambda x: np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)


u1, x, t, time = solver_emv(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
u2, x, t = solver_imsp(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
u3, x, t = solver_cnsp(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
t = T
u_e = np.sin(np.pi*x)*np.exp(-(np.pi**2)*t) \
      + 0.1*np.sin(100*np.pi*x)*np.exp(-(100**2)*(np.pi**2)*t)
# u_i = np.exp(-t)*np.sin(x)
# compare the exact solution with the numerical solution
plt.plot(x, u1, "-.", label="Explicit Euler")
plt.plot(x, u2, "-.", label="Implicit Euler")
plt.plot(x, u3, "-.", label="Crank-Nicolson")
plt.plot(x, u_e, label="Exact solution")
plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.show()
