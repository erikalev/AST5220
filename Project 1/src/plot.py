import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

Mpc = 1e6*3.08567758e16

data = np.loadtxt("data.dat")
dat_int = np.loadtxt("eta_intp.dat")

eta_intp = dat_int

eta = data[:, 0]
x = data[:, 1]
a = np.exp(x)
z = 1./a - 1
Omega_m = data[:, 2]
Omega_b = data[:, 3]
Omega_r = data[:, 4]
Omega_lambda = data[:, 5]
Omega_nu = data[:, 6]
H = data[:, 7]

plt.xlabel("x = ln(a)", size=20)
plt.ylabel("$\eta$(x) [cMpc]", size=20)
plt.semilogy(x, eta/Mpc)
plt.grid()
plt.show()
"""
plt.plot(x, Omega_b, x, Omega_m, x, Omega_r, x, Omega_lambda)
plt.legend(["$\Omega_b$", "$\Omega_m$", "$\Omega_r$", "$\Omega_{\Lambda}$"], loc="center left")
plt.xlabel("x = ln(a)", size=20)
plt.grid()
plt.show()


plt.plot(z, H/1000*Mpc)
ax = plt.gca()
ax.invert_xaxis()
plt.grid()
plt.xlabel("z", size=20)
plt.ylabel("H(z) [km s$^{-1}$ Mpc$^{-1}$]", size=20)
plt.show()

plt.semilogy(x, H/1000*Mpc)
plt.xlabel("x = ln(a)", size=20)
plt.ylabel("H(x) [km s$^{-1}$ Mpc$^{-1}$]", size=20)
plt.grid()
plt.show()
"""

x = np.linspace(x[0], x[-1], len(eta_intp))
plt.semilogy(x, eta_intp/Mpc)
plt.xlabel("x = ln(a)", size=20)
plt.ylabel("$\eta$(x) [cMpc]", size=20)
plt.legend(["interpolated"], loc = "center right")
plt.grid()
plt.show()

