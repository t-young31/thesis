import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

plt.style.use('paper')


def q_t_ew_classical():

    def exp_integrand(x):
        return np.exp(- beta * k * np.exp(a * np.abs(x)))

    integral = integrate.quad(exp_integrand, 0.0, 10.0)[0]
    config_integral = 2.0 * np.exp(beta * k) * integral

    return (np.sqrt(2.0 * m * np.pi / beta) / h_au) * config_integral


if __name__ == '__main__':
    k = 0.00220873
    a = 2.832 * 0.529177
    m = 48 * 1822.888486
    T = 298.15*2

    h_au = 2.0 * np.pi
    kb_au = 3.1668114E-6    # hartrees K-1

    beta = 1.0 / (kb_au * T)   # J -> J mol-1 -> kJ mol-1 -> Ha

    eigenvals = [
        0.000431944,
        0.00103395,
        0.00149685,
        0.00194847,
        0.00236908,
        0.00278874,
        0.00319431,
        0.00360176,
        0.00400121,
        0.00440349,
        0.00480107,
        0.00520147,
        0.00559951,
        0.00600038,
        0.00639984,
        0.00680259,
        0.00720439,
        0.00760902,
        0.00801412,
        0.00842203,
        0.00882995,
        0.00924116
    ]
    ns = list(range(len(eigenvals)+1))

    z_classical = q_t_ew_classical()
    plt.plot(ns, [z_classical for _ in ns])

    z_qunatums = [0.0]
    for eigenval in eigenvals:
        z_qunatum = np.exp(-beta * (eigenval ))
        z_qunatums.append(z_qunatum + z_qunatums[-1])

    plt.plot(ns, z_qunatums)
    plt.savefig('tmp.pdf')
