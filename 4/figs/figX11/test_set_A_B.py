import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
plt.style.use('paper')

j_mol_to_kcal = 1.0 / (1000 * 4.184)


class oneATM:
    plata2015_calc_jkmol = np.array([-171.8, -218.2])

    valsov_calc_jkmol = np.array([-141.7, -141.7, -133.1, -138.1, -142.6,
                                  -99.0, -150.1, -150.2, -108.5, -98.5])


class oneM:

    plata2015_calc_jkmol = np.array([-144.7, -190.8])

    valsov_calc_jkmol = np.array([-113.3, -113.3, -105.6, -111.3, -115.8, -72.2,
                                  -123.4, -125.1, -81.8, -71.9])


class Constants:

    hbar_au = 1.0
    h_au = hbar_au * 2.0 * np.pi
    kb_au = 3.1668114E-6    # hartrees K-1
    h_SI = 6.62607004E-34   # J s

    na = 6.02214086E23       # molecules mol-1
    kb_SI = 1.38064852E-23   # J K-1
    kb_JKmol = kb_SI * na    # J K-1 mol-1
    kb_kcalKmol = kb_SI / (4.184 * 1000)  # kcal K-1 mol-1

    k_b = kb_SI

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

        :param length_au:
        :return:
        """

        q_t = np.sqrt((2.0 * np.pi * self.mass_au * (1.0 / self.beta_au))
                      / (Constants.h_au ** 2)) ** 3 * length_au ** 3
        return Constants.r * (1.5 + np.log(q_t))

    @property
    def s_t_igm_1m(self):
        return self._s_t_igm(length_au=self.l_1molar_au)

    def _q_t_ew(self):
        def exp_integrand(r, beta, a, b):
            return r ** 2 * np.exp(- beta * a * np.exp(b * r))

        cap_lambda = ((2.0 * self.mass_au * np.pi) / (
                    self.beta_au * Constants.h_au ** 2)) ** 1.5
        integral = integrate.quad(exp_integrand, 0.0, 10.0,
                                  args=(self.beta_au, self.a_au, self.b_inv_au))[0]

        return 4.0 * np.pi * np.exp(self.beta_au * self.a_au) * cap_lambda * integral

    @property
    def s_t_ew(self):
        q_t = self._q_t_ew()

        def integrand(r, beta, a, b):
            return r ** 2 * np.exp(-beta * a * (np.exp(b * r) - 1.0) + b * r)

        integral = integrate.quad(integrand, 0.0, 10.0, args=(self.beta_au, self.a_au, self.b_inv_au))[0]

        cap_lambda = ((2.0 * self.mass_au * np.pi) / (self.beta_au * Constants.h_au ** 2)) ** 1.5
        term_4 = 4.0 * np.pi * (self.a_au * self.beta_au * cap_lambda / q_t) * integral

        analytic_s = Constants.r * (1.5 - self.a_au * self.beta_au + np.log(q_t) + term_4)

        return analytic_s

    def sigma_s_t_ew(self, sigma_k_kcal, sigma_a_ang):
        """Standard deviation in the entropy given std. devs. in k and a"""
        eps = 1E-9

        s_curr_k_a = self.s_t_ew
        self.a_au += eps        # Yes k is a..
        ds_dk = (self.s_t_ew - s_curr_k_a) / eps
        self.a_au -= eps

        self.b_inv_au += eps
        ds_da = (self.s_t_ew - s_curr_k_a) / eps
        self.b_inv_au -= eps

        sigma_k = sigma_k_kcal * Constants.kcal_mol_to_au
        sigma_a = sigma_a_ang * Constants.inverse_ang_inverse_au

        return np.sqrt(ds_dk**2 * sigma_k**2 + ds_da**2 + sigma_a**2)

    def __init__(self, mass_amu,
                 k_kcal=0.6183,
                 a_inv_ang=3.10788,
                 temp_K=298.15):
        """

        :param mass_amu:
        :param k_kcal: (average over solute solvent systems methane,
                        acetonitrile etc)
        :param a_inv_ang: (average over solute solvent systems methane,
                         acetonitrile etc)
        :param temp_K:
        """

        self.mass_amu = mass_amu
        self.mass_au = mass_amu * Constants.amu_to_au

        self.a_au = k_kcal * Constants.kcal_mol_to_au
        self.b_inv_au = a_inv_ang * Constants.inverse_ang_inverse_au

        self.temp_K = temp_K
        self.beta_au = 1.0 / (Constants.kb_au * temp_K)

        v_eff_1atm_m3 = Constants.kb_SI * temp_K / Constants.atm_to_pa
        l_1atm = v_eff_1atm_m3 ** (1 / 3) * Constants.m_to_ang
        self.l_1atm_au = l_1atm * Constants.ang_to_au

        v_eff_1molar_m3 = 1.0 / (Constants.n_a * (1.0 / Constants.dm_to_m) ** 3)
        l_1molar = v_eff_1molar_m3 ** (1 / 3) * Constants.m_to_ang
        self.l_1molar_au = l_1molar * Constants.ang_to_au


class Reaction:

    @property
    def dS_trans_1M(self):
        return (sum(prod.s_t_igm_1m for prod in self.prods)
                - sum(reac.s_t_igm_1m for reac in self.reacs))

    @property
    def dS_trans_EW(self):
        return (sum(prod.s_t_ew for prod in self.prods)
                - sum(reac.s_t_ew for reac in self.reacs))

    def dS_EW(self, ds_1m_igm):
        """Subtract the IGM translational entropy and add the EW equivalent"""
        return ds_1m_igm - self.dS_trans_1M + self.dS_trans_EW

    def sigma_dS_EW(self, sigma_k, sigma_a):

        if len(self.prods) != 1:
            raise NotImplementedError

        sigma_s_reacs = np.sqrt(sum(reac.sigma_s_t_ew(sigma_k, sigma_a) ** 2
                                    for reac in self.reacs))

        return np.sqrt(sigma_s_reacs ** 2
                       + self.prods[0].sigma_s_t_ew(sigma_k, sigma_a) ** 2)

    def __init__(self, reacs, prods=None):

        self.reacs = reacs
        if prods is None:
            # Assume the product is just a sum of reactants, with the default
            # a and k
            prods = [Solute(mass_amu=sum(reac.mass_amu for reac in reacs),
                            temp_K=reacs[0].temp_K)]

        self.prods = prods


class EW:

    rmbh1 = Reaction(reacs=[Solute(mass_amu=112.18, temp_K=317.6),
                            Solute(mass_amu=86.09, temp_K=317.6)])

    rmbh2 = Reaction(reacs=[Solute(mass_amu=151.12, temp_K=328.4),
                            Solute(mass_amu=86.09, temp_K=328.4)])

    valsov_temps = [369.15, 369.15, 322.35, 303.15, 303.15, 303.15, 303.15, 252.15, 303.15, 303.15]

    r1 = Reaction(reacs=[Solute(mass_amu=50.49, temp_K=valsov_temps[0]),
                         Solute(mass_amu=58.08, temp_K=valsov_temps[0])])

    # R1 in a different solvent, not just a repeat of the same reaction
    r2 = Reaction(reacs=[Solute(mass_amu=50.49, temp_K=valsov_temps[1]),
                         Solute(mass_amu=58.08, temp_K=valsov_temps[1])])

    r3 = Reaction(reacs=[Solute(mass_amu=141.94, temp_K=valsov_temps[2]),
                         Solute(mass_amu=93.11, temp_K=valsov_temps[2])])

    r4 = Reaction(reacs=[Solute(mass_amu=141.94, temp_K=valsov_temps[3]),
                         Solute(mass_amu=138.10, temp_K=valsov_temps[3])])

    r5 = Reaction(reacs=[Solute(mass_amu=141.94, temp_K=valsov_temps[4]),
                         Solute(mass_amu=155.56, temp_K=valsov_temps[4])])

    r6 = Reaction(reacs=[Solute(mass_amu=141.94, temp_K=valsov_temps[5]),
                         Solute(mass_amu=79.90, temp_K=valsov_temps[5])])

    r7 = Reaction(reacs=[Solute(mass_amu=155.97, temp_K=valsov_temps[6]),
                         Solute(mass_amu=138.10, temp_K=valsov_temps[6])])

    r8 = Reaction(reacs=[Solute(mass_amu=137.02, temp_K=valsov_temps[7]),
                         Solute(mass_amu=109.17, temp_K=valsov_temps[7])])

    r9 = Reaction(reacs=[Solute(mass_amu=256.36, temp_K=valsov_temps[8]),
                         Solute(mass_amu=42.02, temp_K=valsov_temps[8])])

    r10 = Reaction(reacs=[Solute(mass_amu=256.36, temp_K=valsov_temps[9]),
                          Solute(mass_amu=35.45, temp_K=valsov_temps[9])])


def plot_1M():

    plt.scatter(plata_temps*plata2015_expt_jkmol*j_mol_to_kcal,
                plata_temps*oneM.plata2015_calc_jkmol*j_mol_to_kcal,
                label='Set $\\bf{A}$')
    plt.scatter(valsov_temps*valsov_expt_jkmol*j_mol_to_kcal,
                valsov_temps*oneM.valsov_calc_jkmol*j_mol_to_kcal,
                label='Set $\\bf{B}$',
                color='tab:red')

    expt_ts = (plata_temps*j_mol_to_kcal*plata2015_expt_jkmol).tolist() + (valsov_temps*j_mol_to_kcal*valsov_expt_jkmol).tolist()
    cals_ts = (plata_temps*j_mol_to_kcal*oneM.plata2015_calc_jkmol).tolist() + (valsov_temps*j_mol_to_kcal*oneM.valsov_calc_jkmol).tolist()

    print('MAD_IGM = ',
          np.average(np.abs(np.array(expt_ts) - np.array(cals_ts))),
          'kcal mol-1')

    return None


def plot_EW():
    # ------------------------------- Set A ---------------------------------
    EW.plata2015_calc_jkmol = np.array(
        [EW.rmbh1.dS_EW(ds_1m_igm=oneM.plata2015_calc_jkmol[0]),
         EW.rmbh2.dS_EW(ds_1m_igm=oneM.plata2015_calc_jkmol[1])])

    EW.plata2015_calc_std_devs = np.array(
        [EW.rmbh1.sigma_dS_EW(k_sigma, a_sigma),
         EW.rmbh2.sigma_dS_EW(k_sigma, a_sigma)])

    # ------------------------------- Set B ---------------------------------
    EW.valsov_calc_jkmol = np.array(
        [EW.r1.dS_EW(ds_1m_igm=oneM.valsov_calc_jkmol[0]),
         EW.r2.dS_EW(ds_1m_igm=oneM.valsov_calc_jkmol[1]),
         EW.r3.dS_EW(ds_1m_igm=oneM.valsov_calc_jkmol[2]),
         EW.r4.dS_EW(ds_1m_igm=oneM.valsov_calc_jkmol[3]),
         EW.r5.dS_EW(ds_1m_igm=oneM.valsov_calc_jkmol[4]),
         EW.r6.dS_EW(ds_1m_igm=oneM.valsov_calc_jkmol[5]),
         EW.r7.dS_EW(ds_1m_igm=oneM.valsov_calc_jkmol[6]),
         EW.r8.dS_EW(ds_1m_igm=oneM.valsov_calc_jkmol[7]),
         EW.r9.dS_EW(ds_1m_igm=oneM.valsov_calc_jkmol[8]),
         EW.r10.dS_EW(ds_1m_igm=oneM.valsov_calc_jkmol[9])])

    EW.valsov_calc_std_devs = np.array([EW.r1.sigma_dS_EW(k_sigma, a_sigma),
                                        EW.r2.sigma_dS_EW(k_sigma, a_sigma),
                                        EW.r3.sigma_dS_EW(k_sigma, a_sigma),
                                        EW.r4.sigma_dS_EW(k_sigma, a_sigma),
                                        EW.r5.sigma_dS_EW(k_sigma, a_sigma),
                                        EW.r6.sigma_dS_EW(k_sigma, a_sigma),
                                        EW.r7.sigma_dS_EW(k_sigma, a_sigma),
                                        EW.r8.sigma_dS_EW(k_sigma, a_sigma),
                                        EW.r9.sigma_dS_EW(k_sigma, a_sigma),
                                        EW.r10.sigma_dS_EW(k_sigma, a_sigma)])

    plt.errorbar(plata_temps * plata2015_expt_jkmol * j_mol_to_kcal,
                 plata_temps * EW.plata2015_calc_jkmol * j_mol_to_kcal,
                 yerr=plata_temps * EW.plata2015_calc_std_devs * j_mol_to_kcal,
                 linestyle="None",
                 marker='o',
                 label='Set $\\bf{A}$')
    plt.errorbar(valsov_temps * valsov_expt_jkmol * j_mol_to_kcal,
                 valsov_temps * EW.valsov_calc_jkmol * j_mol_to_kcal,
                 yerr=valsov_temps * EW.valsov_calc_std_devs * j_mol_to_kcal,
                 label='Set $\\bf{B}$',
                 linestyle="None",
                 marker='o',
                 color='tab:red')


    expt_ts = (plata_temps*j_mol_to_kcal*plata2015_expt_jkmol).tolist() + (valsov_temps*j_mol_to_kcal*valsov_expt_jkmol).tolist()
    cals_ts = (plata_temps*j_mol_to_kcal*EW.plata2015_calc_jkmol).tolist() + (valsov_temps*j_mol_to_kcal*EW.valsov_calc_jkmol).tolist()

    print('MAD_exp = ',
          np.average(np.abs(np.array(expt_ts) - np.array(cals_ts))),
          'kcal mol-1')


    return None


if __name__ == '__main__':

    ks = [1.048, 0.529, 0.679, 0.545, 0.446, 0.415, 0.53, 1.005, 0.368]
    k_sigma = np.std(np.array(ks))
    as_ = [2.918, 2.793, 2.736, 4.075, 2.930, 3.431, 4.083, 2.127, 2.878]
    a_sigma = np.std(np.array(as_))

    plata_temps = np.array([317.6, 328.4])
    plata2015_expt_jkmol = np.array([-104.2, -129.7])

    valsov_temps = np.array([369.15, 369.15, 322.35, 303.15, 303.15, 303.15,
                             303.15, 252.15, 303.15, 303.15])
    valsov_expt_jkmol = np.array([-44, -69, -24.5, -80.8, -38.7, -28.4, -67.3,
                                  -51, -32.6, -18.4])

    # plot_1M()
    plot_EW()

    min_ts, max_ts = -18, 2

    plt.plot([min_ts, max_ts], [min_ts, max_ts], c='k')
    plt.ylim(min_ts, max_ts)
    plt.xlim(min_ts, max_ts)

    plt.legend()
    plt.xlabel(r'$T\Delta S_{expt}$ / kcal mol$^{-1}$')
    plt.ylabel(r'$T\Delta S_{calc}$ / kcal mol$^{-1}$')

    plt.tight_layout()
    plt.savefig('test_set_A_B_IGM.pdf')
