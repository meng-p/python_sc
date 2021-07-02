from special_topic import solver_em
import matplotlib.pyplot as plt
# import sympy as sym
import numpy as np

Nx = 199
L = 1
T = 0
a = 1
dx = L/Nx
dt = 0.5*(dx**2)
mu = a*dt/(dx**2)  # mu=0.5
# dt = 0.5*(dx**2) # dt=0.005
a = 1
dt = 0.5*(dx**2)
mu = a*dt/(dx**2)  # mu=0.5


def u01(x):
    return np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)
# u01 = lambda x: np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)


u, x, t, time = solver_em(In=u01, alpha=a, L=L, T=T, dx=dx, dt=dt, mu=mu)
plt.plot(x, u, "-", label="t=0")
T = 4.57*10**(-5)
u, x, t, time = solver_em(In=u01, alpha=a, L=L, T=T, dx=dx, dt=dt, mu=mu)
plt.plot(x, u, "-", label="t=4.57*10**(-5)")
T = 2.33*10**(-1)
u, x, t, time = solver_em(In=u01, alpha=a, L=L, T=T, dx=dx, dt=dt, mu=mu)
plt.plot(x, u, "-", label="t=2.33*10**(-1)")
T = 4.67*10**(-1)
u, x, t, time = solver_em(In=u01, alpha=a, L=L, T=T, dx=dx, dt=dt, mu=mu)
plt.plot(x, u, "-", label="t=4.67*10**(-1)")
plt.legend()
plt.xlabel("x")
plt.ylabel("U")
plt.show()
