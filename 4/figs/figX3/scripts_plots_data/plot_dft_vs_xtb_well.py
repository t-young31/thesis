import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
plt.style.use('paper')


def exp_well(x, k, a):
    """V (x) = k(exp(a|x|) − 1)"""
    return k * (np.exp(a*np.abs(x)) - 1.0)


def ho_well(x, k):
    """V (x) = k(exp(a|x|) − 1)"""
    return k * x ** 2 / 2.0


def plot_well(name, color, label):
    """Plot the data plus the fitted well"""

    xs = np.linspace(-1.5, 1.5, num=7)

    ha_to_kcal = 627.509
    all_vs = ha_to_kcal * np.loadtxt(f'{name}.txt')

    vs = np.average(all_vs, axis=0)
    vs_err = np.std(all_vs, axis=0) / np.sqrt(all_vs.shape[0])
    print(np.average(vs_err))

    plt.errorbar(xs, vs,
                 xerr=np.zeros_like(vs_err),
                 yerr=vs_err,
                 fmt='o',
                 color=color,
                 zorder=20)

    more_rs = np.linspace(-1.5, 1.5, num=200)
    # plt.plot(more_rs, exp_well(more_rs, 10, 1.5))

    sigmas = 10 * np.exp(1.5*xs**2)
    # plt.plot(more_rs, 10 * np.exp(0.9 * more_rs**2), c='k', alpha=0.1)

    opt, conv = curve_fit(exp_well, xs, vs,
                          p0=np.array([2, 1.5]),
                          sigma=sigmas
                          )

    plt.plot(more_rs, exp_well(more_rs, *opt), lw=1.4,
             label=label, color=color)

    return vs


if __name__ == '__main__':

    xtb_vs = plot_well(name='GAP_rPBE0-D3_nvt_300K.xyz_small_potentials', color='tab:blue',
              label='GFN2-XTB')
    dft_vs = plot_well(name='GAP_rPBE0-D3_dft_potentials', color='k',
              label='PBE0-D3BJ')

    print(np.mean(np.abs(xtb_vs - dft_vs)))
    exit()

    plt.plot([-2.5, 2.5], [0, 0], ls='-', c='k')

    plt.legend()

    plt.ylabel('V(x) / kcal mol$^{-1}$')
    plt.xlabel('x / Å')
    plt.xlim(-1.6, 1.6)

    plt.ylim(-10, 140)

    plt.savefig('dft_vs_xtb_well.pdf')
