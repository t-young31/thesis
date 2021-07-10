from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('paper')
blues = plt.get_cmap('Blues')


def radial(y, r, m, k, a, l, E):
    psi, phi = y
    return [phi,
            -2.0*phi/r - (2*m*(E - k*(np.exp(a*r) - 1.0)) - (l*(l+1)/r**2))*psi]


if __name__ == '__main__':

    a = 3.0                     # au=1
    k = 0.000054696784          # au
    m = 48 * 1822.888486     # au

    rs = np.linspace(0.0001, 1.7,num=2000)
    sol = odeint(radial,
                 y0=np.array([1.0, -0.00001]),
                 t=rs,
                 args=(m, k, a, 0, 0.000206980766))
    plt.plot(rs, sol[:, 0], lw=2, c=blues(0.5),
             label='$E_0$ for.')

    rs = np.linspace(1.7, 0.0001, num=2000)
    sol = odeint(radial,
                 y0=np.array([0.0001, 0.0]),
                 t=rs,
                 args=(m, k, a, 0, 0.000207029908))
    plt.plot(rs, sol[:, 0]/sol[-20, 0], lw=2,
             c=blues(0.9),
             ls='--',
             label='$E_0$ rev.')

    plt.ylim(-0.25, 1.05)
    # plt.xlim(0, 0.1)
    plt.legend()
    plt.ylabel(r'$\tilde{\varphi}_{l=0}(r)$')
    plt.xlabel('$r$ / au')
    plt.tight_layout()
    plt.savefig('first_eigenfunction.pdf')
