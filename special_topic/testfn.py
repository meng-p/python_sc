from special_topic import solver_emv
import matplotlib.pyplot as plt
# import sympy as sym
import numpy as np
# from explicit import solver_em

dx = 0.01
L = 1
T = 0.1
Nx = int(round(L/dx))
# dt = 0.5*(dx**2) # dt=0.005
a = 1
dt = 0.5*(dx**2)
mu = a*dt/(dx**2)  # mu=0.5
x = np.linspace(0, L, Nx+1)
u0 = np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)
# plt.plot(x, u0, '-', label='inital')


def u01(x):
    return np.sin(np.pi*x) + 0.1*np.sin(100*np.pi*x)


u, x, t, cpu = solver_emv(In=u01, alpha=a, L=L, T=T, dx=dx, dt=dt, mu=mu)

t = T
u_e = np.sin(np.pi*x)*np.exp(-(np.pi**2)*t) \
      + 0.1*np.sin(100*np.pi*x)*np.exp(-(100**2)*(np.pi**2)*t)

print(cpu)
plt.plot(x, u0, label="Intial condition")
plt.legend()
plt.show()
plt.plot(x, u, "-.", label="Finite difference approximation")
plt.plot(x, u_e, label="Exact solution")
plt.legend()
plt.show()
