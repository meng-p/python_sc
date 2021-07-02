from special_topic import solver_em2
import matplotlib.pyplot as plt
# import sympy as sym
import numpy as np
# from mpl_toolkits import mplot3d

fig = plt.figure()
x = np.linspace(0, np.pi, 20)  # x \in [0,pi], Nx+1=20
t = np.linspace(0, 3, 40)  # t \in [0,3], Nt+1=40
x, t = np.meshgrid(x, t)

Z = np.exp(-t)*np.sin(x)
ax = plt.axes(projection='3d')
ax.plot_surface(x, t, Z, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')  # 3d plot
ax.set_xlabel("x")

ax.set_ylabel("t")

ax.set_zlabel("u")
plt.show()


def u01(x):  # initial condition
    return np.sin(x)


L = np.pi
T = 3
Nx = 14  # choose a different step size
Nt = 199
a = 1
dx = L/Nx
dt = T/Nt
mu = dt/(dx**2)
u, x, t = solver_em2(In=u01, alpha=a, L=L, T=T, dx=dx, dt=dt, mu=mu)

x, t = np.meshgrid(x, t)
ax = plt.axes(projection='3d')
ax.plot_surface(x, t, u, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
# ax.set_title('surface')
ax.set_xlabel("x")

ax.set_ylabel("t")

ax.set_zlabel("U")
plt.show()

Nx = 14
Nt = 94
a = 1
dx = L/Nx
dt = T/Nt
mu = dt/(dx**2)
u, x, t = solver_em2(In=u01, alpha=a, L=L, T=T, dx=dx, dt=dt, mu=mu)

x, t = np.meshgrid(x, t)
ax = plt.axes(projection='3d')
ax.plot_surface(x, t, u, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_xlabel("x")

ax.set_ylabel("t")

ax.set_zlabel("U")
plt.show()
