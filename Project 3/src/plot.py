import numpy as np
import matplotlib.pyplot as plt
import itertools

time_data = np.loadtxt("data_time.dat")

eta = time_data[:, 0]
x1 = time_data[:, 1]
Omega_m = time_data[:, 2]
Omega_b = time_data[:, 3]
Omega_r = time_data[:, 4]
Omega_lambda = time_data[:, 5]
Omega_nu = time_data[:, 6]
#plt.plot(x1, Omega_m, x1, Omega_r)
#plt.show()

data        = np.loadtxt("data.dat")
x       = np.loadtxt("x_data.dat")
n = 1500
k = 6

#x         = np.zeros(n)
Phi         = np.zeros((n, k))
dPhi        = np.zeros((n, k))
Psi         = np.zeros((n, k))
dPsi        = np.zeros((n, k))
delta       = np.zeros((n, k))
ddelta      = np.zeros((n, k))
delta_b     = np.zeros((n, k))
ddelta_b    = np.zeros((n, k))
v           = np.zeros((n, k))
dv          = np.zeros((n, k))
v_b         = np.zeros((n, k))
dv_b        = np.zeros((n, k))
Theta0      = np.zeros((n, k))
dTheta0     = np.zeros((n, k))
Theta1      = np.zeros((n, k))
dTheta1     = np.zeros((n, k))
Theta2      = np.zeros((n, k))
dTheta2     = np.zeros((n, k))
Theta3      = np.zeros((n, k))
dTheta3     = np.zeros((n, k))
Theta4      = np.zeros((n, k))
dTheta4     = np.zeros((n, k))

n_value = 0
k_value = 0

for i in range(n*k):
    Phi[n_value, k_value]         = data[i, 0]
    dPhi[n_value, k_value]        = data[i, 1]
    Psi[n_value, k_value]         = data[i, 2]
    dPsi[n_value, k_value]        = data[i, 3]
    delta[n_value, k_value]       = data[i, 4]
    ddelta[n_value, k_value]      = data[i, 5]
    delta_b[n_value, k_value]     = data[i, 6]
    ddelta_b[n_value, k_value]    = data[i, 7]
    v[n_value, k_value]           = data[i, 8]
    dv[n_value, k_value]          = data[i, 9]
    v_b[n_value, k_value]         = data[i, 10]
    dv_b[n_value, k_value]        = data[i, 11]
    Theta0[n_value, k_value]      = data[i, 12]
    dTheta0[n_value, k_value]     = data[i, 13]
    Theta1[n_value, k_value]      = data[i, 14]
    dTheta1[n_value, k_value]     = data[i, 15]
    Theta2[n_value, k_value]      = data[i, 16]
    dTheta2[n_value, k_value]     = data[i, 17]
    Theta3[n_value, k_value]      = data[i, 18]
    dTheta3[n_value, k_value]     = data[i, 19]
    Theta4[n_value, k_value]      = data[i, 20]
    dTheta4[n_value, k_value]     = data[i, 21]
    k_value +=1
    if k_value == k:
        k_value = 0
        n_value += 1

x_eq = np.zeros(len(x))
x_eq[:] = -7.9
y_eq = np.zeros(len(x))
y_eq = np.linspace(1.5, max(delta[-1, :]), len(x))




x_tick = [-20, -15, -10, -7.9, -5, 0]
for i in range(k):
    plt.semilogy(x, delta[:, i])
plt.legend(["kc/H$_0$ = 0.1", "kc/H$_0$ = 8.36", "kc/H$_0$ = 85.90", "kc/H$_0$ = 245.1", "kc/H$_0$ = 636.8", "kc/H$_0$ = 1000.0"], loc = "upper left")
plt.semilogy(x_eq, y_eq, dashes=[6, 2])
plt.grid()
plt.xticks(x_tick)
plt.ylabel("$\\delta$", size=20)
ax = plt.gca()
ax.set_xticks(x_tick)
ax.set_xticklabels([-20, -15, -10, "x$_{eq}$", -5, 0])
plt.xlabel("x", size=20)
plt.show()

#plt.subplot(122)
for i in range(k):
    plt.plot(x, delta_b[:, i])
plt.legend(["kc/H$_0$ = 0.1", "kc/H$_0$ = 8.36", "kc/H$_0$ = 85.90", "kc/H$_0$ = 245.1", "kc/H$_0$ = 636.8", "kc/H$_0$ = 1000.0"], loc = "upper left")
plt.grid()
plt.ylim(1e-3, 1e5)
plt.yscale("log")
plt.xlabel("x", size=20)
plt.ylabel("$\\delta_b$", size=20)
plt.show()
"""
for i in range(k):
    plt.plot(x, v[:, i])
plt.legend(["kc/H$_0$ = 0.1", "kc/H$_0$ = 8.36", "kc/H$_0$ = 85.90", "kc/H$_0$ = 245.1", "kc/H$_0$ = 636.8", "kc/H$_0$ = 1000.0"], loc = "upper left")
plt.grid()
plt.xlabel("x", size=20)
plt.ylabel("v", size=20)
plt.show()

for i in range(k):
    plt.plot(x, v_b[:, i])
plt.legend(["kc/H$_0$ = 0.1", "kc/H$_0$ = 8.36", "kc/H$_0$ = 85.90", "kc/H$_0$ = 245.1", "kc/H$_0$ = 636.8", "kc/H$_0$ = 1000.0"], loc = "upper left")
plt.grid()
plt.xlabel("x", size=20)
plt.ylabel("v$_b$", size=20)
plt.show()

for i in range(k):
    plt.plot(x, Phi[:, i])
plt.legend(["kc/H$_0$ = 0.1", "kc/H$_0$ = 8.36", "kc/H$_0$ = 85.90", "kc/H$_0$ = 245.1", "kc/H$_0$ = 636.8", "kc/H$_0$ = 1000.0"], loc = "upper left")
plt.grid()
plt.xlabel("x", size=20)
plt.ylabel("$\\Phi$", size=20)
plt.show()

for i in range(k):
    plt.plot(x, Psi[:, i])
plt.grid()
plt.ylim(-1.0, 0.2)
plt.legend(["kc/H$_0$ = 0.1", "kc/H$_0$ = 8.36", "kc/H$_0$ = 85.90", "kc/H$_0$ = 245.1", "kc/H$_0$ = 636.8", "kc/H$_0$ = 1000.0"], loc = "upper left")
plt.xlabel("x", size=20)
plt.ylabel("$\\Psi$", size=20)
plt.show()


for i in range(k):
    plt.plot(x, Theta0[:, i])
plt.grid()
plt.ylim([-1.5, 1.5])
plt.legend(["kc/H$_0$ = 0.1", "kc/H$_0$ = 8.36", "kc/H$_0$ = 85.90", "kc/H$_0$ = 245.1", "kc/H$_0$ = 636.8", "kc/H$_0$ = 1000.0"], loc = "upper left")
plt.xlabel("x", size=20)
plt.ylabel("$\\Theta_0$", size=20)
plt.show()

for i in range(k):
    plt.plot(x, Theta1[:, i])
    plt.ylim([-1.5, 1.5])
plt.legend(["kc/H$_0$ = 0.1", "kc/H$_0$ = 8.36", "kc/H$_0$ = 85.90", "kc/H$_0$ = 245.1", "kc/H$_0$ = 636.8", "kc/H$_0$ = 1000.0"], loc = "upper left")
plt.grid()
plt.xlabel("x", size=20)
plt.ylabel("$\\Theta_1$", size=20)
plt.show()

for i in range(k):
    plt.plot(x, Theta2[:, i])
    plt.ylim([-1.5, 1.5])
plt.legend(["kc/H$_0$ = 0.1", "kc/H$_0$ = 8.36", "kc/H$_0$ = 85.90", "kc/H$_0$ = 245.1", "kc/H$_0$ = 636.8", "kc/H$_0$ = 1000.0"], loc = "upper left")
plt.grid()
plt.xlabel("x", size=20)
plt.ylim([-0.5, 0.5])
plt.ylabel("$\\Theta_2$", size=20)
plt.show()

for i in range(k):
    plt.plot(x, Theta3[:, i])
    plt.ylim([-1.5, 1.5])
plt.legend(["kc/H$_0$ = 0.1", "kc/H$_0$ = 8.36", "kc/H$_0$ = 85.90", "kc/H$_0$ = 245.1", "kc/H$_0$ = 636.8", "kc/H$_0$ = 1000.0"], loc = "upper left")
plt.grid()
plt.xlabel("x", size=20)
plt.ylabel("$\\Theta_3$", size=20)
plt.show()
"""
