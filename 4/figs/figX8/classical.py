import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.stats import linregress
plt.style.use('paper')


def integrand(r, _beta, _k, _a):
    return r**2 * np.exp(- _beta * _k * np.exp(_a * r))


def z_classical():

    z_c = (np.exp(beta * k)
           * (np.sqrt((2.0 * m * np.pi * kb_T)) / h)**3
           * 4 * np.pi * quad(integrand, 0.0, 10.0,
                              args=(beta, k, a))[0])

    return z_c


def plot_integrand():

    rs = np.linspace(0, 2, num=500)

    fig, ax = plt.subplots()
    ax_2 = ax.twinx()

    ax.plot(rs, integrand(rs, beta, k, a) * 10 ** 3, label='integral', lw=1.5)
    ax.set_xlabel('$r$ / au')
    ax.set_ylabel('integral / au${}^3 \\times 10^{-3}$')
    ax.legend(loc=(0.75, 0.1))

    ax_2.plot(rs, k * (np.exp(a * rs) - 1.0) / kb_T,
              label='V(r)',
              c='k',
              lw=2)
    ax_2.set_ylim(-2.5, 50)
    ax_2.set_ylabel('$V(r)$ / $k_B T$')
    ax_2.set_xticks([0, 0.5, 1.0, 1.5, 2.0])

    ax_2.legend(loc=(0.75, 0.2))

    plt.tight_layout()
    plt.savefig('classical_integrand.pdf')

    exit()
    return None


if __name__ == '__main__':

    k = 0.00220873             # Ha
    m = 48 * 1822.888486       # m
    a = 2.832 * 0.529177  # a0^-1
    kb = 3.1668114E-6     # Ha K-1
    T = 298
    h = 2.0 * np.pi       # Ha s

    kb_T = kb * T
    beta = 1.0 / (kb * T)

    # plot_integrand()

    rs = np.linspace(0, 2, num=500)

    fig, ax_k = plt.subplots(figsize=(5.5, 5))
    ax_a = ax_k.twiny()

    ks = np.linspace(0.001, 0.003, num=10)
    integrands_k = [quad(integrand, 0.0, 5.0, args=(beta, k, a))[0]
                    for k in ks]

    ax_k.scatter(ks/kb_T, np.log10(np.array(integrands_k)),
                 c='tab:blue',
                 label='k')
    ax_k.legend(loc=(0.85, 0.9))
    ax_k.set_xlabel('$k$ / $k_B T$')
    ax_k.set_yticks([-2, -2.5, -3, -3.5])
    # ax_k.set_xticks([1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3])

    lin = linregress(ks/kb_T, np.log10(np.array(integrands_k)))
    more_ks = np.linspace(np.min(ks/kb_T), np.max(ks/kb_T), num=200)
    ax_k.plot(more_ks, lin.slope * more_ks + lin.intercept,
              ls='--')

    as_ = np.linspace(2.0 * 0.529177, 5.0 * 0.529177, num=10)
    integrands_a = [quad(integrand, 0.0, 2.0, args=(beta, k, a))[0]
                    for a in as_]

    ax_a.scatter(as_, np.log10(np.array(integrands_a)),
                 c='tab:orange',
                 label='a')
    lin = linregress(as_, np.log10(np.array(integrands_a)))
    more_as = np.linspace(np.min(as_), np.max(as_), num=200)
    ax_a.plot(more_as, lin.slope * more_as + lin.intercept,
              ls='--',
              c='tab:orange')

    ax_a.set_xlabel('$a$ / au$^{-1}$')
    ax_a.legend(loc=(0.85, 0.82))

    ax_k.set_ylabel('log$_{10}$(integral / au$^3$)')

    plt.tight_layout()
    plt.savefig('integrand_a_k.pdf')
