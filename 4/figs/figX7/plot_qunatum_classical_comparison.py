import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
plt.style.use('paper')


def integrand(r):
    return r**2 * np.exp(- beta * k * np.exp(a * r))


def z_classical():

    # plt.plot(rs, integrand(rs))
    # plt.plot(rs, k * (np.exp(a * rs) - 1.0))
    # plt.savefig('tmp_int.pdf')

    z_c = (np.exp(beta * k)
           * (np.sqrt((2.0 * m * np.pi * kb_T)) / h)**3
           * 4 * np.pi * quad(integrand, 0.0, 10.0)[0])

    return z_c


if __name__ == '__main__':

    limiting = True

    if limiting:                   # Classical limit of the problem
        k = 0.0001                 # Ha
        m = 48 * 1822.888486 * 10  # m_e

    else:                          # CO2
        k = 0.00220873             # Ha
        m = 48 * 1822.888486       # m_e

    a = 2.832 * 0.529177  # a0^-1

    kb = 3.1668114E-6     # Ha K-1
    T = 50
    h = 2.0 * np.pi       # Ha s

    kb_T = kb * T
    beta = 1.0 / (kb * T)

    z_quantum = [0]

    for line in open('eigenvals.txt', 'r').readlines()[1:]:
        eigenval, l_num = line.split('\t')

        e_j = (2*float(l_num) + 1) * np.exp(-(float(eigenval))/kb_T)
        z_quantum.append(z_quantum[-1] + e_j)

    n_eigenvals = len(z_quantum)
    plt.plot(list(range(n_eigenvals)), z_quantum,
             label='$Z_q$')
    plt.plot(list(range(n_eigenvals)), n_eigenvals * [z_classical()],
             label='$Z_c$')
    plt.xlabel('$\\tilde{n}$')
    plt.ylabel('$Z$')
    plt.legend()
    plt.tight_layout()
    plt.savefig('quantum_classical_comparison_limiting.pdf')

    print(z_classical() / z_quantum[-1])
