from special_topic import solver_em, solver_cn1, solver_im
# import sympy as sym
import matplotlib.pyplot as plt
import numpy as np

dx = 0.1  # set the step size in space
L = 1
T = 0.1  # a given final time
dt = 0.5*(dx**2)  # the step size in time depends on dx
mu = dt/(dx**2)  # mu=0.5


def u01(x):  # initial condition
    return np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)
# u01 = lambda x: np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)


u1, x, t = solver_im(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
u2, x, t, time = solver_em(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
u3, x, t = solver_cn1(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
t = T
u_i = np.sin(np.pi*x)*np.exp(-(np.pi**2)*t) \
      + 0.1*np.sin(100*np.pi*x)*np.exp(-(100**2)*(np.pi**2)*t)

# compare the exact solution with the numerical solution
plt.plot(x, u2, "-.", label="Explicit Euler")
plt.plot(x, u1, "-.", label="Implicit Euler")
plt.plot(x, u3, "-.", label="Crank-Nicolson")
plt.plot(x, u_i, label="Exact solution")
plt.legend()
plt.show()
