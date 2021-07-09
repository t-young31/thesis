import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
plt.style.use('paper')
blues = plt.get_cmap('Blues')


def gto0(x):
    return 0.19682158E-01 * np.exp(-13.0107010 * x**2)


def gto1(x):
    return 0.13796524 * np.exp(-1.9622572 * x**2)


def gto2(x):
    return 0.47831935 * np.exp(-0.44453796 * x**2)


def s_gtos(x):
    return gto0(x) + gto1(x) + gto2(x)


def s_true(x):
    return 2 * np.exp(-np.abs(x))


if __name__ == '__main__':

    # N_true = np.sqrt(1.0 / quad(lambda x: s_true(x)**2 * x**2, 0, np.inf)[0])
    N_true = 1.0

    # N_gtos = np.sqrt(1 / quad(lambda x: s_gtos(x)**2 * x**2, 0, np.inf)[0])
    N_gtos = 2.6929627898296626

    xs = np.linspace(-6, 6, num=500)
    plt.plot(xs, N_true * s_true(xs), lw=1.5, c='k', label='exp(-r)')


    plt.plot(xs, N_gtos * gto0(xs), lw=1.0, c=blues(0.3),
             # label='$c_1$ exp(-α$_1r$)'
             )

    plt.plot(xs, N_gtos * gto1(xs), lw=1.0, c=blues(0.3),
             # label='$c_2$ exp(-α$_2r$)'
             )

    plt.plot(xs, N_gtos * gto2(xs), lw=1.0, c=blues(0.3),
             # label='$c_3$ exp(-α$_3r$)'
             )

    plt.plot(xs, N_gtos * s_gtos(xs), lw=1.5, c=blues(0.8),
             label='$\Sigma_{i=1}^3 c_i$ exp(-α$_ir$)')

    plt.ylim(0, 2.2)
    plt.xlim(-6, 6)
    plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0])
    plt.legend()
    plt.ylabel('$R_{1s}(x)$')
    plt.xlabel(r'$x\; /\; a_0$')
    plt.tight_layout()
    plt.savefig('gto_AO.pdf')
