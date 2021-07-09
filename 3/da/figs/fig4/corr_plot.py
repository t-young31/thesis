import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
plt.style.use('paper2')


# quinone: 1  12  14  15  16  21  G3  G4  G7  G8
dE_bind_C1 = [-4.0, -5.3, -7.5, -6.8, -8.4,  8.1, -9.5, -14.6, -8.6, -5.0]
dE_bind_C1_expt = [-5.3, None, -7.4, None, None, None, -10.5, -12.2, -3.6, -3.6]
dE_bind_C2 = [-3.3, -5.3, -6.1, -6.6, -6.4, 5.1, -6.8, -11.2, -7.1, -3.4]
dE_bind_C2_expt = [-4.1, -5.0, -4.7, -6.3, -5.7, 0.0, -5.4, -8.3, None, None]


def plot_de_correlation():

    all_expt = dE_bind_C1_expt + dE_bind_C2_expt
    all_calc = dE_bind_C1 + dE_bind_C2

    all_expt_not_none = [val for val in all_expt if val is not None]
    all_calc_not_none = []
    for i in range(len(all_expt)):
        if all_expt[i] is not None:
            all_calc_not_none.append(all_calc[i])

    min_e, max_e = -16, 10
    m, c, r, p, err = linregress(all_expt_not_none, all_calc_not_none)
    print('MSE = ', np.round(np.average(np.array(all_calc_not_none) - np.array(all_expt_not_none)), 2),
          'MAD = ', np.round(np.average(np.abs(np.array(all_calc_not_none) - np.array(all_expt_not_none))), 2),
          'r^2 =', np.round(np.square(r), 3))

    plt.scatter(dE_bind_C2_expt, dE_bind_C2, marker='o', s=80, edgecolors='k', linewidth='0.5',label='C-2')
    plt.scatter(dE_bind_C1_expt, dE_bind_C1, marker='o', s=80, edgecolors='k', linewidth='0.5',label='C-1')
    plt.xlabel('expt. $\Delta G_{bind}$ / kcal mol$^{-1}$')
    plt.ylabel('calcd. $\Delta E_{bind}$ / kcal mol$^{-1}$')
    plt.ylim(min_e, max_e)
    plt.xlim(min_e, max_e)
    # plt.plot([min_e, max_e], [min_e * m + c, max_e * m + c], c='k', ls='--')
    # plt.fill_between([min_e, max_e], [min_e - 1, max_e - 1], [min_e + 1, max_e + 1], color='k', alpha=0.1)
    plt.plot([min_e, max_e], [min_e - 1, max_e - 1], color='k', lw=0.3)

    x_range = np.array([min_e, max_e])
    plt.plot(x_range, m * x_range + c, color='k', lw=1, ls='--')


    plt.plot([min_e, max_e], [min_e + 1, max_e + 1], color='k', lw=0.3)
    plt.tight_layout()

    return plt.show()


def plot_dde_correlation():

    min_e, max_e = -16, 10

    dde_bind_C2_expt = [val - dE_bind_C1_expt[0] if val is not None else None for val in dE_bind_C2_expt]
    dde_bind_C2 = [val - dE_bind_C1[0] for val in dE_bind_C2]
    dde_bind_C1_expt = [val - dE_bind_C1_expt[0] if val is not None else None for val in dE_bind_C1_expt]
    dde_bind_C1 = [val - dE_bind_C1[0] for val in dE_bind_C1]

    all_expt = dde_bind_C1_expt + dde_bind_C2_expt
    all_calc = dde_bind_C1 + dde_bind_C2

    all_expt_not_none = [val for val in all_expt if val is not None]
    all_calc_not_none = []
    for i in range(len(all_expt)):
        if all_expt[i] is not None:
            all_calc_not_none.append(all_calc[i])

    m, c, r, p, err = linregress(all_expt_not_none, all_calc_not_none)
    print('MSE = ', np.round(np.average(np.array(all_calc_not_none) - np.array(all_expt_not_none)), 2),
          'MAD = ', np.round(np.average(np.abs(np.array(all_calc_not_none) - np.array(all_expt_not_none))), 2),
          'r^2 =', np.round(np.square(r), 3))

    plt.scatter(dde_bind_C2_expt, dde_bind_C2, marker='o', s=80, edgecolors='k', linewidth='0.5', label='C-2')
    plt.scatter(dde_bind_C1_expt, dde_bind_C1, marker='o', s=80, edgecolors='k', linewidth='0.5', label='C-1')
    plt.xlabel('expt. $\Delta\Delta G_{bind}$ / kcal mol$^{-1}$')
    plt.ylabel('calcd. $\Delta\Delta E_{bind}$ / kcal mol$^{-1}$')
    plt.ylim(min_e, max_e)
    plt.xlim(min_e, max_e)

    x_range = np.array([min_e, max_e])
    plt.plot(x_range, m * x_range + c, color='k', lw=1, ls='--')
    # plt.plot([min_e, max_e], [min_e * m + c, max_e * m + c], c='k', ls='--')
    # plt.fill_between([min_e, max_e], [min_e - 1, max_e - 1], [min_e + 1, max_e + 1], color='k', alpha=0.1)
    plt.plot([min_e, max_e], [min_e - 1, max_e - 1], color='k', lw=0.3)
    # plt.plot([min_e, max_e], [min_e, max_e], color='k', lw=0.3)
    plt.plot([min_e, max_e], [min_e + 1, max_e + 1], color='k', lw=0.3)
    plt.tight_layout()

    return plt.show()


if __name__ == '__main__':

    plot_de_correlation()
    plot_dde_correlation()
