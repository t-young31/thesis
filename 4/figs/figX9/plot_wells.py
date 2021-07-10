import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
plt.style.use('paper')
ha_to_kcal = 627.509


def exp_well(x, k, a, c):
    """V (x) = k(exp(a|x|) − 1)"""
    return k * (np.exp(a*np.abs(x)) - 1.0) + c


def plot_well(name, ax_):
    """Plot the data plus the fitted well"""
    all_vs = ha_to_kcal * np.loadtxt(f'{name}.txt')

    if len(xs) != all_vs.shape[1]:
        raise ValueError(name + f" {len(xs)} != {all_vs.shape[1]}")

    vs = np.average(all_vs, axis=0)

    vs_err = np.std(all_vs, axis=0) / np.sqrt(all_vs.shape[0])

    ax_.errorbar(xs, vs,
                 xerr=np.zeros_like(vs_err),
                 yerr=vs_err,
                 fmt='o',
                 color='k',
                 zorder=20)

    more_rs = np.linspace(np.min(xs), np.max(xs), num=200)
    # plt.plot(more_rs, exp_well(more_rs, 10, 1.5))

    # sigmas = 10 * np.exp(3.5*xs**2)
    # plt.plot(more_rs, 10 * np.exp(0.9 * more_rs**2), c='k', alpha=0.1)

    opt, conv = curve_fit(exp_well, xs[vs < 10], vs[vs < 10],
                          p0=np.array([2, 1.5, 0.0]),
                          # sigma=sigmas
                          )

    print("-in-".join(name.split('_')[:2]),
          f'k = {np.round(opt[0], 3)} kcal mol-1',
          f'a = {np.round(opt[1], 3)} Å-1')

    ax_.plot(more_rs, exp_well(more_rs, *opt), lw=1.4)
    ax_.set_xlim(-1.1, 1.,1)
    ax_.set_ylim(-1, 10)
    if 'water' in name:
        ax_.set_ylabel('V(x) / kcal mol$^{-1}$')
    if 'alanine' in name:
        ax_.set_xlabel('x / Å')
    ax_.plot([-3, 3], [0, 0], ls='-', c='k')

    return None


if __name__ == '__main__':

    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(10.5, 9),
                           sharex='col',
                           sharey='row')

    xs = np.linspace(-1, 1, num=10)         # water

    plot_well(name='methane_water',
              ax_=ax[0, 0])

    plot_well(name='methane_mecn',
              ax_=ax[0, 1])

    plot_well(name='methane_benzene',
              ax_=ax[0, 2])

    plot_well(name='co2_water',
              ax_=ax[1, 0])

    plot_well(name='co2_mecn',
              ax_=ax[1, 1])

    plot_well(name='co2_benzene',
              ax_=ax[1, 2])

    plot_well(name='alanine_water',
              ax_=ax[2, 0])

    plot_well(name='alanine_mecn',
              ax_=ax[2, 1])

    plot_well(name='alanine_benzene',
              ax_=ax[2, 2])

    plt.tight_layout()
    plt.savefig('wells.pdf')
