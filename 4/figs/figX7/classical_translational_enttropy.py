from constants import *
from scipy import integrate
import matplotlib.pyplot as plt


def s_infinite_spherical_well(temp_K, mass_au, c_au):

    z = ((2.0 * np.pi * mass_au * Constants.kb_au * temp_K) / (Constants.h_au ** 2))**1.5 * ((4.0 * np.pi * c_au**3) / 3.0)

    return 1.5 * Constants.kb_JKmol + Constants.kb_JKmol * np.log(z)


def s_pib(temp_K, mass_au, l_au):

    z = ((2.0 * np.pi * mass_au * Constants.kb_au * temp_K) / (Constants.h_au**2))**1.5 * l_au**3

    return 1.5 * Constants.kb_JKmol + Constants.kb_JKmol * np.log(z)


def exp_well_potential(r_au, k_au, a_au):
    return k_au * np.exp(a_au * r_au)


def s_exp_well(temp_K, mass_au, a_au, k_au):

    def config_integral(r, beta_au, k_au, a_au):
        return r**2 * np.exp(- beta_au * exp_well_potential(r, k_au, a_au))

    beta_au = 1.0 / (Constants.kb_au * temp_K)

    integral = integrate.quad(config_integral, 0.0, 20.0, args=(beta_au, k_au, a_au))[0]

    z = (4.0 * np.pi * ((2.0 * np.pi * mass_au * Constants.kb_au * temp_K) / (Constants.h_au**2))**1.5 *
         np.exp(beta_au * k_au) *
         integral)

    print(z)
    return 1.5 * Constants.kb_JKmol + Constants.kb_JKmol * np.log(z)


def plot_potentials(c_au, k_au, a_au):

    xs = np.linspace(0, 6, 200)
    plt.plot(xs, [exp_well_potential(x, k_au, a_au) for x in xs])
    plt.plot([0, c_au-1E-8, c_au], [0.0, 0.0, 1E4])
    plt.ylim(-1, 4)

    return plt.show()


if __name__ == '__main__':

    # c = 2.0 * 1.88                      # 2 Å to Bohr
    # a = 3.0                             # au
    # k = 0.000054696784                  # au
    temp = 298.15                       # K
    m = 48 * Constants.amu_to_au     # amu to au

    # plot_potentials(c, k, a)

    v_eff_1atm_m3 = Constants.k_b * temp / Constants.atm_to_pa
    l_1atm_au = v_eff_1atm_m3**(1/3) * Constants.m_to_bohr

    v_eff_1molar_m3 = 1.0 / (Constants.n_a * (1.0 / Constants.dm_to_m)**3)
    l_1molar_au = v_eff_1molar_m3**(1/3) * Constants.m_to_bohr

    l = l_1atm_au
    # m = 44.01 * Constants.amu_to_au     # CO2
    # m = 89.09 * Constants.amu_to_au     # alanine
    # m = 78.11 * Constants.amu_to_au     # benzene

    # k_kjmol, a_ang = 1.186, 4.92          # CO2 in H2O
    # k_kjmol, a_ang = 0.467, 5.5343        # CO2 in MeOH
    # k_kjmol, a_ang = 2.7, 2.2             # CO2 in DCM

    # k_kjmol, a_ang = 20, 3.62             # alanine in H2O
    # k_kjmol, a_ang = 11.56, 3.60          # alanine in Meoh
    # k_kjmol, a_ang = 16.8, 4.67           # alanine in DCM

    # k_kjmol, a_ang = 1.68, 5.01             # benzene in H2O
    # k_kjmol, a_ang = 1.41, 4.94          # benzene in Meoh
    # k_kjmol, a_ang = 7.0, 2.01           # benzene in DCM


    # k, a = k_kjmol * Constants.kj_mol_to_au, a_ang * Constants.inverse_ang_inverse_au

    k, a = 0.0001, 2.832 * 0.529177

    print('S_PIB       = ', s_pib(temp, m, l_au=l))
    print('S_spherical = ', s_infinite_spherical_well(temp, m, c_au=l/2.0))
    print('S_exp       = ', s_exp_well(temp, m, a_au=a, k_au=k))
