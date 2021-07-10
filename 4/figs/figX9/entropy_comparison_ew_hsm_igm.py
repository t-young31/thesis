"""
Calculate the translational entropy with EW, HSM and IGM models

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
plt.style.use("paper")
ha2kjmol = 627.5 * 4.184


class Constants:

    hbar_au = 1.0
    h_au = hbar_au * 2.0 * np.pi
    kb_au = 3.1668114E-6    # hartrees K-1
    h_SI = 6.62607004E-34   # J s

    na = 6.02214086E23       # molecules mol-1
    kb_SI = 1.38064852E-23   # J K-1
    kb_JKmol = kb_SI * na    # J K-1 mol-1
    kb_kcalKmol = kb_SI / (4.184 * 1000)  # kcal K-1 mol-1

    k_b = kb_kcalKmol

    n_a = 6.022140857E23                        # molecules mol-1
    r = k_b * n_a                               # J K-1 mol-1

    h = 6.62607004E-34                          # J s
    atm_to_pa = 101325                          # Pa
    dm_to_m = 0.1                               # m
    amu_to_kg = 1.660539040E-27                 # Kg
    c = 299792458                               # m s-1
    c_in_cm = c * 100                           # cm s-1
    ang_to_m = 1E-10                            # m
    ang_to_au = 1.88973                         # au Å-1
    m_to_ang = 1E10                             # Å
    m_to_bohr = 1.89e+10                        # au m-1
    amu_to_au = 1822.888486                     # m_e amu-1
    kj_mol_to_au = 0.00038087980                # Ha (kJ mol-1)-1
    kcal_mol_to_au = 0.001593601                # Ha (kcal mol-1)-1
    inverse_ang_inverse_au = 1.0 / 1.88973      # au-1 Å


class Solute:

    def _s_t_igm(self, length_au):
        """
        S = k_B T dln(q_t)/dT + k_B ln(q_t)
          = n k_B ((T / q) (3q / 2T) + ln(q))                     dq/dT = 3q / 2T
          = n k_B (3/2 + ln(q))

        :param mass_au:
        :param temp_K:
        :param length_au:
        :return:
        """

        q_t = np.sqrt((2.0 * np.pi * self.mass_au * Constants.kb_au * temp_K) / (
                    Constants.h_au ** 2)) ** 3 * length_au ** 3
        return Constants.r * (1.5 + np.log(q_t))

    @property
    def s_t_igm_1atm(self):
        return self._s_t_igm(length_au=self.l_1atm_au)

    @property
    def s_t_igm_1m(self):
        return self._s_t_igm(length_au=self.l_1molar_au)

    def s_t_hsm(self, omega_au):
        """
        S = k_B T dln(q_t)/dT + k_B ln(q_t)
          = (n k_B T) d ln(q)/dT + n k_B ln(q)
          = n k_B (3 hbar omega beta coth(hbar omega beta / 2) / 2T + ln(q))

        :param temp_K:
        :param omega_au:
        :return:
        """
        raise NotImplementedError

        beta_au = 1.0 / (Constants.kb_au * temp_K)
        q_t = 1.0 / (2.0 * np.sinh(
            Constants.hbar_au * omega_au * beta_au / 2.0)) ** 3
        term1 = 3.0 * (Constants.hbar_au * omega_au * beta_au / 2.0) / np.tanh(
            Constants.hbar_au * omega_au * beta_au / 2.0)

        return Constants.r * (term1 + np.log(q_t))

    def _q_t_ew(self):
        def exp_integrand(r, beta, a, b):
            return r ** 2 * np.exp(- beta * a * np.exp(b * r))

        cap_lambda = ((2.0 * self.mass_au * np.pi) / (
                    beta_au * Constants.h_au ** 2)) ** 1.5
        integral = integrate.quad(exp_integrand, 0.0, 10.0,
                                  args=(beta_au, self.a_au, self.b_inv_au))[0]

        return 4.0 * np.pi * np.exp(beta_au * self.a_au) * cap_lambda * integral

    @property
    def s_t_ew(self):
        beta_au = 1.0 / (Constants.kb_au * temp_K)
        q_t = self._q_t_ew()

        def integrand(r, beta, a, b):
            return r ** 2 * np.exp(-beta * a * (np.exp(b * r) - 1.0) + b * r)

        integral = integrate.quad(integrand, 0.0, 10.0, args=(beta_au, self.a_au, self.b_inv_au))[0]

        cap_lambda = ((2.0 * self.mass_au * np.pi) / (beta_au * Constants.h_au ** 2)) ** 1.5
        term_4 = 4.0 * np.pi * (self.a_au * beta_au * cap_lambda / q_t) * integral

        analytic_s = Constants.r * (1.5 - self.a_au * beta_au + np.log(q_t) + term_4)

        # d_temp = 1E-10
        # dlnq_dtemp = ((np.log(
        #     q_t_ew(1.0 / (Constants.kb_au * (temp_K + d_temp)), self.a_au, self.b_inv_au,
        #            mass_au)) - np.log(q_t))
        #               / d_temp)

        # numerical_s = Constants.r * (temp_K * dlnq_dtemp + np.log(q_t))
        # print('Numerical derivative / analytic derivative = ', numerical_s / analytic_s)

        return analytic_s

    def __init__(self, mass_amu, k_kcal, a_inv_ang):

        self.mass_au = mass_amu * Constants.amu_to_au

        self.a_au = k_kcal * Constants.kcal_mol_to_au
        self.b_inv_au = a_inv_ang * Constants.inverse_ang_inverse_au

        # Harmonic oscillator
        # k_au = k_kjmol * Constants.kj_mol_to_au
        # omega_au = np.sqrt(k_au / mass_au)

        v_eff_1atm_m3 = Constants.kb_SI * temp_K / Constants.atm_to_pa
        l_1atm = v_eff_1atm_m3 ** (1 / 3) * Constants.m_to_ang
        self.l_1atm_au = l_1atm * Constants.ang_to_au

        v_eff_1molar_m3 = 1.0 / (Constants.n_a * (1.0 / Constants.dm_to_m) ** 3)
        l_1molar = v_eff_1molar_m3 ** (1 / 3) * Constants.m_to_ang
        self.l_1molar_au = l_1molar * Constants.ang_to_au


def ST(S):
    return np.round(temp_K * S, decimals=3)


if __name__ == '__main__':

    temp_K = 298.15
    beta_au = 1.0 / (Constants.kb_au * temp_K)

    Methane_Water = Solute(mass_amu=16.04, k_kcal=1.048, a_inv_ang=2.918)
    Methane_Acetonitrile = Solute(mass_amu=16.04, k_kcal=0.529, a_inv_ang=2.793)
    Methane_Benzene = Solute(mass_amu=16.04, k_kcal=0.679, a_inv_ang=2.736)

    CO2_Water = Solute(mass_amu=44.01, k_kcal=0.545, a_inv_ang=4.075)
    CO2_Acetonitrile = Solute(mass_amu=44.01, k_kcal=0.446 , a_inv_ang=2.93)
    CO2_Benzene = Solute(mass_amu=44.01, k_kcal=0.415, a_inv_ang=3.431)

    Alanine_Water = Solute(mass_amu=89.09, k_kcal=0.53, a_inv_ang=4.083)
    Alanine_Acetonitrile = Solute(mass_amu=89.09, k_kcal=1.005, a_inv_ang=2.127)
    Alanine_Benzene = Solute(mass_amu=89.09, k_kcal=0.368, a_inv_ang=2.878)

    systems = [Methane_Water, Methane_Acetonitrile, Methane_Benzene,
               CO2_Water, CO2_Acetonitrile, CO2_Benzene,
               Alanine_Water, Alanine_Acetonitrile, Alanine_Benzene]

    rels = []
    for system in systems:
        rel = ST(system.s_t_ew) / ST(system.s_t_igm_1m)
        rels.append(rel)

        print(ST(system.s_t_igm_1atm),
              ST(system.s_t_igm_1m),
              ST(system.s_t_ew),

              'kcal mol-1',
              sep='  &  ')


    print(np.average(np.array(rels)),
          '±', np.std(np.array(rels))/np.sqrt(len(rels)))
