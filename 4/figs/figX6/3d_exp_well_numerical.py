from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np


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


def radial_h_potential(y, r, m, l, E):
    psi, phi = y
    return [phi,  -2.0*phi/r - (2*m*(E + 1/r) - (l*(l+1)/r**2))*psi]


def h_atom():
    """H atom solved numerically"""

    m = 1
    eigval = -1/2

    rs = np.linspace(0.000001, 5.0, num=2000)
    sol = odeint(radial_h_potential,
                 y0=np.array([1, -0.001]),
                 t=rs,
                 args=(m, 0, eigval))
    plt.plot(rs, sol[:, 0])
    plt.plot(rs, np.exp(-rs))
    plt.savefig('tmp.pdf')
    return


if __name__ == '__main__':

    # h_atom()

    a = 3.0                     # au=1
    k = 0.000054696784          # au
    # temp = 298.15               # K
    m = 48 * 1822.888486     # au

    rs = np.linspace(1.7, 0.0001, num=2000)
    # rs = np.linspace(0.0001, 1.7, num=2000)

    range_1 = [0.000207029879, 0.000207030293]
    range_2 = [0.000477839, 0.000477870]
    range_3 = [0.000781328436, 0.000781359045]
    range_4 = [0.001115069891, 0.001115120906]

    range_default = (0.0008, 0.0013)

    # E_1 = 0.000206980766  Ha  from forwards
    # E_1 = 0.000207029908  Ha from backwards

    for E in np.linspace(*range_4, num=100):
        sol = odeint(radial,
                     y0=np.array([0.0001, 0.0]),
                     t=rs,
                     args=(m, k, a, 0, E))
        print(f'E = {E:.12f}    ψ(r=0) =', sol[-1][0])

        # plt.plot(rs[:-20], sol[:-20, 0]/sol[-20, 0])
        plt.plot(rs, sol[:, 0])

    # plt.ylim(-1.1, 1.1)
    # plt.xlim(0, 0.1)
    plt.savefig('tmp.pdf')
