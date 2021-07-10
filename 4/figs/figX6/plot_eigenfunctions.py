from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
plt.style.use('paper')
reds = plt.get_cmap('Reds')

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=[reds(i/5+0.1) for i in range(4)])


def radial(y, r, m, k, a, l, E):
    """
    Radial ODE with an exponential well

    φ  = dψ/dr

    dφ/dr = -2φ/r - (2m(E-ke^ar) - (l(l+1)/r^2))ψ

    :param y: (np.ndarray) First order ODEs
    --------------------- all in atomic units ------------------
    :param r: (float) Distance
    :param m: (float) Mass
    :param k: (float) Exponential well constant
    :param a: (float) Exponential well exponent
    :param E: (float) Energy

    :return: (np.ndarray): Step of the ODE
    """
    psi, phi = y

    dydt = [phi,   # φ  = dψ/dr
            -2.0*phi/r - (2*m*(E - k*(np.exp(a*r) - 1.0)) - (l*(l+1)/r**2))*psi
            ]
    return dydt


if __name__ == '__main__':

    # h_atom()

    a = 3.0                     # au=1
    k = 0.000054696784          # au
    # temp = 298.15               # K
    m = 48 * 1822.888486     # au

    rs = np.linspace(1.7, 0.0001, num=2000)
    # rs = np.linspace(0.0001, 1.7, num=2000)

    E_0 = (0.000207029879 + 0.000207030293) / 2.0
    E_1 = (0.000477839 + 0.000477870) / 2.0
    E_2 = (0.000781328436 + 0.000781359045) / 2.0
    E_3 = (0.001115069891 + 0.001115120906) / 2.0

    for i, E in enumerate([E_0, E_1, E_2, E_3]):
        print(E)
        sol = odeint(radial,
                     y0=np.array([0.0001, 0.0]),
                     t=rs,
                     args=(m, k, a, 0, E))
        # Normalise to the not unstable first point
        plt.plot(rs[:-20], sol[:-20, 0]/sol[-20, 0],
                 lw=2,
                 label=f'$E_{i}$')

    # plt.ylim(-1.1, 1.1)
    # plt.xlim(0, 0.1)
    plt.legend()
    plt.ylabel(r'$\tilde{\varphi}_{l=0}(r)$')
    plt.xlabel('$r$ / au')
    plt.tight_layout()
    plt.savefig('eigenfunctions.pdf')
