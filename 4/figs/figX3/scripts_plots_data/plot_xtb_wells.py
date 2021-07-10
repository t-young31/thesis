import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
plt.style.use('paper')


def exp_well(x, k, a):
    """V (x) = k(exp(a|x|) − 1)"""
    return k * (np.exp(a*np.abs(x)) - 1.0)


def plot_well(name, label, color, max_fit_energy=500):
    """Plot the data plus the fitted well"""

    xs = np.linspace(-2, 2, num=10)

    ha_to_kcal = 627.509
    all_vs = ha_to_kcal * np.loadtxt(f'{name}.txt')

    vs = np.average(all_vs, axis=0)
    vs_err = np.std(all_vs, axis=0) / np.sqrt(all_vs.shape[0])

    plt.errorbar(xs, vs,
                 xerr=np.zeros_like(vs_err),
                 yerr=vs_err,
                 label=label,
                 fmt='o',
                 color=color)

    more_rs = np.linspace(-2, 2, num=200)
    # plt.plot(more_rs, exp_well(more_rs, 10, 1.5))

    sigmas = 10 * np.exp(1.5*xs**2)
    # plt.plot(more_rs, 10 * np.exp(0.9 * more_rs**2), c='k', alpha=0.1)

    opt, conv = curve_fit(exp_well,
                          xdata=[x for i, x in enumerate(xs)
                                 if vs[i] < max_fit_energy],
                          ydata=[v for i, v in enumerate(vs)
                                 if vs[i] < max_fit_energy],
                          p0=np.array([2, 1.5]),
                          sigma=[s for i, s in enumerate(sigmas)
                                 if vs[i] < max_fit_energy]
                          )

    print(label,
          f'k = {np.round(opt[0], 3)} kcal mol-1',
          f'a = {np.round(opt[1], 3)} Å-1',
          sep='\t')

    plt.plot(more_rs, exp_well(more_rs, *opt),lw=1.2, color=color)

    return None


if __name__ == '__main__':

    plot_well(name='dftb_300_nvt.xyz_potentials', label='DFTB',
              color='tab:blue')
    plot_well(name='GAP_rPBE0-D3_nvt_300K.xyz_potentials', label='revPBE0-D3',
              color='k')
    plot_well(name='tip4p_300_nvt.xyz_potentials', label='MM-TIP4P',
              color='tab:green')

    plt.plot([-2.5, 2.5], [0, 0], ls='-', c='k')

    plt.ylabel('$V(x)$ / kcal mol$^{-1}$')
    plt.xlabel('$x$ / Å')
    plt.legend()
    plt.xlim(-2.1, 2.1)

    plt.ylim(-20, 620)
    plt.tight_layout()
    plt.savefig('xtb_wells.pdf')

    plt.ylim(-10, 100)
    plt.ylabel(None)
    plt.tight_layout()
    plt.savefig('xtb_wells_zoom.pdf')
