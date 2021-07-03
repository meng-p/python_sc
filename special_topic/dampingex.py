from special_topic import solver_em
import matplotlib.pyplot as plt
# import sympy as sym
import numpy as np

# consider the explicit Euler scheme
Nx = 199
L = 1
# first we consider the numerical solution at time t=0
# (matches the initial condition)
T = 0
dx = L/Nx  # step size in space
dt = 0.5*(dx**2)  # step size in time (related to dx)
mu = dt/(dx**2)  # mu=0.5
# dt = 0.5*(dx**2) # dt=0.005


def u01(x):
    return np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)
# u01 = lambda x: np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)


u, x, t, time = solver_em(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
plt.plot(x, u, "-", label="t=0")
# then plot the numerical solution at time T=4.57*10^(-5)
T = 4.57*10**(-5)
u, x, t, time = solver_em(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
plt.plot(x, u, "-", label="t=4.57*10**(-5)")
# plot the numerical solution at time T=2.33*10^(-1)
T = 2.33*10**(-1)
u, x, t, time = solver_em(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
plt.plot(x, u, "-", label="t=2.33*10**(-1)")
# plot the numerical solution at time T=4.67*10^(-1)
T = 4.67*10**(-1)
u, x, t, time = solver_em(In=u01, L=L, T=T, dx=dx, dt=dt, mu=mu)
plt.plot(x, u, "-", label="t=4.67*10**(-1)")
plt.legend()
plt.xlabel("x")
plt.ylabel("U")
plt.show()
