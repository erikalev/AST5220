import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

Mpc = 1e6*3.08567758e16

def plot_project_1():
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
    

    x = np.linspace(x[0], x[-1], len(eta_intp))
    plt.semilogy(x, eta_intp/Mpc)
    plt.xlabel("x = ln(a)", size=20)
    plt.ylabel("$\eta$(x) [cMpc]", size=20)
    plt.legend(["interpolated"], loc = "center right")
    plt.grid()
    plt.show()

def plot_project_2():
    data = np.loadtxt("data.dat")
    x_rec = data[:, 0]
    X_e = data[:, 1]

    tau = data[:, 2]
    tau2 = abs(data[:, 3])
    tau22 = abs(data[:, 4])
    g = data[:, 5]
    dg_dx = data[:, 6]
    dg2_d2x = data[:, 7]
    a = np.exp(x_rec)
    z = 1./a - 1
    
    ax = plt.gca()
    ax.invert_xaxis()
    my_yticks = [10**(-3), 10**(-1), 10**(-1), 1]
    plt.ylim([10**(-4), 1.1])
    plt.xticks([1800, 1400, 1000, 600, 200], size=20)
    plt.xlim([1800, 0])
    plt.yticks([1, 0.1, 0.01, 0.001], size=20)
    plt.xlabel("z", size=20)
    plt.ylabel("X$_e$", size=20)
    plt.semilogy(z, X_e)
    plt.grid()
    plt.show()
    
    plt.xlabel("x = $\ln(a)$", size=20)
    plt.xticks(np.arange(-17.5, 0, 2.5), size=20)
    plt.xlim([-17.5, 0])
    plt.yticks(size=20)
    plt.ylabel("$\\tau$, $|\\tau '|$, $|\\tau ''|$", size=20)
    plt.semilogy(x_rec, tau, x_rec, tau2, x_rec, tau22)
    plt.legend(["$\\tau$", "$|\\tau '|$", "$|\\tau ''|$"])
    plt.grid()
    plt.show()
    
    plt.plot(x_rec, g, x_rec, dg_dx/10.0, x_rec, dg2_d2x/300.0)
    plt.xlabel("x = $\ln(a)$", size=20)
    plt.ylabel("$\\tilde{g}$, $\\tilde{g}'$, $\\tilde{g}$''", size=20)
    plt.xticks([-7.4, -7.2, -7.0, -6.8, -6.6, -6.4, -6.2, -6], size=20)
    plt.xlim([-7.5, -6])
    plt.ylim([-6, 6])
    plt.yticks([-2, 0, 2, 4])
    plt.legend(["$\\tilde{g}$", "$\\tilde{g}'$", "$\\tilde{g}$''"])
    plt.grid()
    plt.show()
    
plot_project_1()
#plot_project_2()
