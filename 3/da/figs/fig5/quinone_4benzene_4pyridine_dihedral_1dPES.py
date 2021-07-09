"""
Data from N_Pd_C_C_scan.out, which performed an unrelaxed scan of a

----------------------------------------------------------
                  benzene
                   /
                  /
benzene  ---  quinone  ---- benzene
                /
               /
            benzene
----------------------------------------------------------
and
----------------------------------------------------------

                  pyridine
                   /
                  /
pyridine  ---  quinone  ---- pyridine
                /
               /
            pyridine
----------------------------------------------------------
systems

Calculated at M06-2X/def2-TZVP. Charges are calculated with the HIRSHFELD scheme and correspond to the sum of the
quinone atoms
"""


import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
plt.style.use('paper')


dihedral_deg = [-45, -41.9, -38.79, -35.69, -32.59, -29.48, -26.38, -23.28, -20.17, -17.07, -13.97, -10.86, -7.76,
                -4.66, -1.55, 1.55, 4.66, 7.76, 10.86, 13.97, 17.07, 20.17, 23.28, 26.38, 29.48, 32.59, 35.69, 38.79,
                41.9, 45]


pyridine_rel_energies_kcal = [1.705266497, 1.710972229, 1.685555868, 1.615677311, 1.473753051, 1.287334385, 1.083012855,
                              0.864434412, 0.669605896, 0.499578308, 0.367638608, 0.253493879, 0.151024585, 0.060677844,
                              0, 2.17679E-05, 0.060528349, 0.151311497, 0.253740776, 0.367683161, 0.499539208,
                              0.669596026, 0.864252248, 1.083085858, 1.287333939, 1.47372866, 1.615688073, 1.685574016,
                              1.710859963, 1.705200603]  # C-C dist in cage = 10.256 Å actual 10.252

pyridine_charges = [-0.02464, -0.024285, -0.02368, -0.023399, -0.024002, -0.025421, -0.027653, -0.030152, -0.032761,
                    -0.036029, -0.039269, -0.04206, -0.044129, -0.045515, -0.046195, -0.046195, -0.045517, -0.044131,
                    -0.042061, -0.03927, -0.036028, -0.032763, -0.030153, -0.027653, -0.025421, -0.024001, -0.023398,
                    -0.02368, -0.024285, -0.024641]  # Total quinone partial atomic charge

pyridine_n_charges = [-0.75427, -0.754146, -0.753892, -0.753326, -0.752817, -0.751985, -0.750969, -0.749823, -0.748627,
                      -0.747644, -0.746888, -0.746183, -0.745467, -0.745141, -0.745003, -0.745003, -0.74514, -0.745467,
                      -0.746185, -0.746888, -0.747645, -0.748625, -0.749823, -0.750969, -0.751985, -0.752817, -0.753327,
                      -0.753891, -0.754148, -0.75427]  # Sum of nitrogen charges

onering_dihedral_deg = [0.0, 6.207, 12.414, 18.621, 24.828, 31.035, 37.242, 43.449, 49.656, 55.863, 62.07, 68.277,
                        74.484, 80.691, 86.898, 93.105]

onepyridine_rel_energies = [-2.302512224, -2.221333163, -2.133312069, -1.979151399, -1.692701285, -1.404858793,
                            -1.175709212, -1.050975744, -0.884237154, -0.76582919, -0.747494814, -0.741281544,
                            -0.723122834, -0.640412623, -0.642987398, -0.642911905]
#                          Relative energy zeroed to benzene + pyridine + quinone at infinite separation

benzene_rel_energies_kcal = [0, 0.018793568, 0.053649549, 0.099668498, 0.137708891, 0.223423352, 0.337083217,
                             0.434982542, 0.553956711, 0.691966997, 0.840896698, 0.99644239, 1.139933512, 1.228321354,
                             1.249913102, 1.248605078, 1.228610607, 1.140860536, 0.996646083, 0.841228383, 0.692286771,
                             0.554902028, 0.434208954, 0.337882891, 0.220951786, 0.137393309, 0.100326689, 0.054834426,
                             0.018530684, 0.001167388]  # N-N dist in cage = 10.405Å, actual 10.393

benzene_charges = [0.018701, 0.019144, 0.020382, 0.02187, 0.023245, 0.024046, 0.023554, 0.022088, 0.020937, 0.01939,
                   0.017743, 0.015899, 0.014528, 0.013667, 0.013171, 0.01317, 0.013663, 0.014528, 0.015898, 0.017744,
                   0.019388, 0.020939, 0.022093, 0.023554, 0.024044, 0.023246, 0.021869, 0.020384, 0.019141, 0.018703]
#                  Total quinone partial atomic charge

benzene_h_charges = [-0.368241, -0.368059, -0.367801, -0.367637, -0.367053, -0.366177, -0.365203, -0.364402, -0.364058,
                     -0.363892, -0.363805, -0.363439, -0.36297, -0.362549, -0.362284, -0.362285, -0.362549, -0.36297,
                     -0.363442, -0.363806, -0.363895, -0.364058, -0.364401, -0.3652, -0.366176, -0.367051, -0.367638,
                     -0.367799, -0.368061, -0.368239]
#                   Sum of all benzene H charges (where H is directed toward the centre

onebenzene_rel_energies = [0.175953572, 0.159512313, -0.015621801, -0.183284545, -0.239555556, -0.295851938,
                           -0.31890044, -0.357568132, -0.307727768, -0.272233325, -0.290989306, -0.310271294,
                           -0.313206951, -0.257078469, -0.272485594, -0.269888751]
#                          Relative energy zeroed to benzene + pyridine + quinone at infinite separation


def cos_func(x, a, b, c, d):
    x_rad = np.deg2rad(x)
    return a * np.cos(b * x_rad + d) + c


def poly_fit(some_x, some_y, order=10):

    more_x = np.linspace(min(some_x), max(some_x), 200)
    poly_fit = np.polyfit(some_x, some_y, order)
    y = np.poly1d(poly_fit)(more_x)

    return more_x, y


def plot_1dpes(units='kJ'):

    if units == 'kcal':
        pyridine_rel_energies = np.array(pyridine_rel_energies_kcal)
        benzene_rel_energies = np.array(benzene_rel_energies_kcal)

    if units == 'kJ':
        pyridine_rel_energies = 4.184 * np.array(pyridine_rel_energies_kcal)
        benzene_rel_energies = 4.184 * np.array(benzene_rel_energies_kcal)


    more_dihedral = np.linspace(min(dihedral_deg), max(dihedral_deg), 200)
    pyridine_cos_params = curve_fit(cos_func, dihedral_deg, pyridine_rel_energies, p0=[-0.7, 5, 0.7, 0.0])[0]
    benzene_cos_params = curve_fit(cos_func, dihedral_deg, benzene_rel_energies, p0=[0.7, 5, -0.7, 0.0])[0]

    plt.scatter(dihedral_deg, pyridine_rel_energies, marker='+', s=60, label='N')
    plt.plot(more_dihedral, [cos_func(x, *pyridine_cos_params) for x in more_dihedral], ls='-', lw=1.5)

    plt.scatter(dihedral_deg, benzene_rel_energies, marker='+', s=60, label='CH')
    plt.plot(more_dihedral, [cos_func(x, *benzene_cos_params) for x in more_dihedral], ls='-', lw=1.5)

    plt.xlabel('$\phi$ / deg')
    plt.ylabel('$\Delta E$ / ' + units + ' mol$^{-1}$')
    plt.tight_layout()
    # plt.legend(prop={'size': 13})

    return plt.show()


def plot_quinone_charges():

    rel_pyridine_charges = np.array(pyridine_charges) - pyridine_charges[0]
    rel_benzene_chagres = np.array(benzene_charges) - benzene_charges[0]

    plt.scatter(dihedral_deg, rel_pyridine_charges, marker='+', label='N')
    plt.plot(*poly_fit(dihedral_deg, rel_pyridine_charges, order=10), ls='--')

    plt.scatter(dihedral_deg, rel_benzene_chagres, marker='+', label='CH')
    plt.plot(*poly_fit(dihedral_deg, rel_benzene_chagres, order=10), ls='--')

    plt.xlabel('$\phi$ / deg')
    plt.ylabel('$\Delta q_{tot}$ / $e$')
    plt.ylim(-0.027, 0.012)
    plt.tight_layout()
    # plt.legend()
    return plt.show()


def plot_x_charges():

    rel_pyridine_charges = np.array(pyridine_n_charges) - pyridine_n_charges[0]
    rel_benzene_chagres = np.array(benzene_h_charges) - benzene_h_charges[0]

    plt.scatter(dihedral_deg, rel_pyridine_charges, marker='+', label='N')
    plt.scatter(dihedral_deg, rel_benzene_chagres, marker='+', label='CH')
    plt.xlabel('$\\theta$ / deg')
    plt.ylabel('$\Delta q$ / $e$')
    plt.legend()
    return plt.show()


def tanh_func(x, a, b, c, d):
    x_rad = np.deg2rad(x)
    return a * np.tanh(b * x_rad + d) + c


def plot_1dpes_1ring():

    pyridine_rel_energies_kJmol = 4.184 * np.array(onepyridine_rel_energies)
    benzene_rel_energies_kJmol = 4.184 * np.array(onebenzene_rel_energies)

    more_dihedral = np.linspace(min(onering_dihedral_deg), max(onering_dihedral_deg), 200)

    pyridine_cos_params = curve_fit(tanh_func, onering_dihedral_deg, pyridine_rel_energies_kJmol, p0=[-0.7, 5, 0.7, 0.0])[0]
    benzene_cos_params = curve_fit(tanh_func, onering_dihedral_deg, benzene_rel_energies_kJmol, p0=[-1, 10, -0.1, -2.0])[0]

    plt.plot([-1, 95], [0, 0], c='k', lw=0.5)

    plt.scatter(onering_dihedral_deg, pyridine_rel_energies_kJmol, marker='+', label='N')
    plt.plot(more_dihedral, [tanh_func(x, *pyridine_cos_params) for x in more_dihedral], ls='--')

    plt.scatter(onering_dihedral_deg, benzene_rel_energies_kJmol, marker='+', label='CH')
    plt.plot(more_dihedral, [tanh_func(x, *benzene_cos_params) for x in more_dihedral], ls='--')

    plt.xlim(min(onering_dihedral_deg), max(onering_dihedral_deg))
    plt.xlabel('$\phi$ / deg')
    plt.ylabel('$\Delta E$ / kJ mol$^{-1}$')
    # plt.legend()

    return plt.show()


if __name__ == '__main__':

    # plot_1dpes(units='kcal')
    plot_quinone_charges()
    # plot_x_charges()
    # plot_1dpes_1ring()
