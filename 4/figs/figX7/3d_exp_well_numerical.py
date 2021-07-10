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


def h_atom(n=2):
    """H atom solved numerically

    Works for 1s, 2s, 2p
    """

    m = 1
    eigval = -1/(2 * n**2)

    rs = np.linspace(0.01, 10.0, num=2000)
    sol = odeint(radial_h_potential,
                 y0=np.array([0.015, 1E-6]),
                 t=rs,
                 args=(m, 1, eigval))

    plt.plot(rs, sol[:, 0])                 # 2nd column is dφ/dr
    plt.plot(rs, np.exp(-rs), label='1s')
    plt.plot(rs, (2 - rs)*np.exp(-rs/2)/2, label='2s')
    plt.plot(rs, rs*np.exp(-rs/2), label='2p')

    plt.legend()
    plt.savefig('tmp.pdf')
    exit()
    return


def old_params():
    #     n l  respectively, indexed from 0 with
    """
    a = 3.0                     # au=1
    k = 0.000054696784          # au
    # temp = 298.15               # K
    m = 48 * 1822.888486     # au
    """
    range_00 = [0.000207029879, 0.000207030293]
    range_10 = [0.000477839, 0.000477870]
    range_20 = [0.000781328436, 0.000781359045]
    range_30 = [0.001115069891, 0.001115120906]

    range_01 = [0.000337342910, 0.000337342922]
    range_11 = [0.000623858533, 0.000623858545]
    range_21 = [0.000942477455, 0.000942477467]

    return


if __name__ == '__main__':

    # h_atom()

    k = 0.00220873
    a = 2.832 * 0.529177
    m = 48 * 1822.888486

    # rs = np.linspace(1.7, 0.003, num=2000)
    rs = np.linspace(0.001, 1.7, num=2000)

    range_default = (0.001537, 0.001538)

    # E_1 = 0.000206980766  Ha  from forwards
    # E_1 = 0.000207029908  Ha from backwards

    for E in np.linspace(*range_default, num=100):
        sol = odeint(radial,
                     y0=np.array([1, 0.0]),
                     t=rs,
                     args=(m, k, a, 1, E))
        print(f'E = {E:.12f}    ψ(r=0) =', sol[-1][0])

        # plt.plot(rs[:-20], sol[:-20, 0]/sol[-20, 0])
        plt.plot(rs, sol[:, 0])

    plt.ylim(-1.1, 1.1)
    # plt.xlim(0, 0.1)
    plt.savefig('tmp.pdf')
