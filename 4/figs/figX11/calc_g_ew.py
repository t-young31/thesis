from otherm import Molecule
import matplotlib.pyplot as plt
from scipy import integrate
import numpy as np
plt.style.use('paper')
reds = plt.get_cmap('Reds')
blues = plt.get_cmap('Blues')

temp = 298.15


class Constants:

    hbar_au = 1.0
    h_au = hbar_au * 2.0 * np.pi
    kb_au = 3.1668114E-6    # hartrees K-1
    h_SI = 6.62607004E-34   # J s

    na = 6.02214086E23       # molecules mol-1
    kb_SI = 1.38064852E-23   # J K-1
    kb_kcalKmol = kb_SI / (4.184 * 1000)  # kcal K-1 mol-1

    j_mol_to_kcal = 1.0 / (4.184 * 1000)    # kcal mol-1 (J mol-1)^-1
    au_to_kcal = 627.509

    n_a = 6.022140857E23                        # molecules mol-1

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
    au_to_j_mol = 1000.0 / kj_mol_to_au         # J mol-1 Ha^-1

    inverse_ang_inverse_au = 1.0 / 1.88973      # au-1 Å

    r = kb_SI * n_a                             # J K mol-1


class EWMoleucle(Molecule):

    @property
    def mass_au(self):
        return (self.mass / Constants.amu_to_kg) * Constants.amu_to_au

    @property
    def s_t_igm(self):
        """
        S = k_B T dln(q_t)/dT + k_B ln(q_t)
          = n k_B ((T / q) (3q / 2T) + ln(q))                  dq/dT = 3q / 2T
          = n k_B (3/2 + ln(q))

        :return:
        """
        effective_volume = (Constants.m_to_bohr**3
                            / (Constants.n_a * (1.0 / Constants.dm_to_m)**3))

        q_t = np.sqrt((2.0 * np.pi * self.mass_au * (1.0 / self.beta_au))
                      / (Constants.h_au ** 2)) ** 3 * effective_volume
        return Constants.r * (2.5 + np.log(q_t))

    def _q_t_ew(self):
        def exp_integrand(r, beta, a, b):
            return r ** 2 * np.exp(- beta * a * np.exp(b * r))

        cap_lambda = ((2.0 * self.mass_au * np.pi) / (
                self.beta_au * Constants.h_au ** 2)) ** 1.5
        integral = integrate.quad(exp_integrand, 0.0, 10.0,
                                  args=(
                                  self.beta_au, self.k_au, self.a_inv_au))[0]

        return 4.0 * np.pi * np.exp(
            self.beta_au * self.k_au) * cap_lambda * integral

    @property
    def s_t_ew(self):
        q_t = self._q_t_ew()

        def integrand(r, beta, a, b):
            return r ** 2 * np.exp(-beta * a * (np.exp(b * r) - 1.0) + b * r)

        integral = integrate.quad(integrand, 0.0, 10.0,
                                  args=(self.beta_au,
                                        self.k_au,
                                        self.a_inv_au))[0]

        cap_lambda = ((2.0 * self.mass_au * np.pi)
                      / (self.beta_au * Constants.h_au ** 2)) ** 1.5
        term_4 = 4.0 * np.pi * (self.k_au * self.beta_au * cap_lambda / q_t) * integral

        analytic_s = Constants.r * (
                    1.5 - self.k_au * self.beta_au + np.log(q_t) + term_4)

        return analytic_s

    @property
    def internal_e_ew(self):
        """
        Calculate the internal energy from an exponeital well rather than a
        standard PIB treatment

        :return: (float)
        """
        cap_lambda = ((2.0 * self.mass_au * np.pi)
                      / (self.beta_au * Constants.h_au ** 2)) ** 1.5

        def integrand(r, beta, a, b):
            return r ** 2 * np.exp(-beta * a * (np.exp(b * r) - 1.0) + b * r)

        integral = integrate.quad(integrand, 0.0, 10.0,
                                  args=(self.beta_au,
                                        self.k_au,
                                        self.a_inv_au))[0]

        term3 = Constants.au_to_j_mol * ((4.0 * np.pi * self.k_au * cap_lambda
                                          / self._q_t_ew())
                                         * integral)

        return 1.5 * Constants.r - self.k_au * Constants.au_to_j_mol + term3

    def calculate_thermochemistry(self, method='grimme', shift=100, w0=100,
                                  alpha=4, calc_sym=True, symm_n=None):

        # Call the super class function for the standard IGM (qRRHO) things
        super().calculate_thermochemistry(temp, '1M', method, shift, w0, alpha,
                                          calc_sym, symm_n)
        return

    def __init__(self, filename,
                 k_kcal=0.6183,
                 a_inv_ang=3.10788):

        super().__init__(filename, is_ts=False, real_freqs=False)

        self.filename = filename
        self.k_au = k_kcal * Constants.kcal_mol_to_au
        self.a_inv_au = a_inv_ang * Constants.inverse_ang_inverse_au

        self.beta_au = 1.0 / (Constants.kb_au * temp)


if __name__ == '__main__':

    monomer = EWMoleucle(filename='monomer.out')
    monomer.calculate_thermochemistry()

    monomer_ew = EWMoleucle(filename='monomer.out')
    monomer_ew.calculate_thermochemistry(method='no_s_trans')


    octamer = EWMoleucle(filename='octamer.out')
    octamer.calculate_thermochemistry()

    octamer_ew = EWMoleucle(filename='octamer.out')

    # Delete the three lowest energy modes from the EW monomer
    for _ in range(3):
        del octamer_ew.freqs[-7]

    octamer_ew.calculate_thermochemistry(method='no_s_trans')

    dE_opt_to_sp = [-0.010322268710041271, -0.010912749090920215, -0.01371703320600659, -0.01564209978593567, -0.01605363347297839, -0.012402004561948843, -0.01375865564099854,-0.01152977633691421]
    dE_opt_to_sp = Constants.au_to_kcal * np.array(dE_opt_to_sp)

    dgs = []
    for i in range(1, 9):
        septamer = EWMoleucle(filename=f'septamer{i}.out')
        septamer.calculate_thermochemistry()

        TdS_trans = Constants.j_mol_to_kcal * temp * (septamer.s_t_ew - octamer_ew.s_t_ew)
        dU_trans = Constants.j_mol_to_kcal*(septamer.internal_e_ew - octamer_ew.internal_e_ew)
        print(dU_trans)
        print(Constants.j_mol_to_kcal*temp*monomer.s_rot )
        exit()

        # TdS_trans = Constants.j_mol_to_kcal * temp * (septamer.s_t_ew + monomer.s_t_ew - octamer_ew.s_t_ew)
        # dU_trans = Constants.j_mol_to_kcal*(septamer.internal_e_ew + monomer.internal_e_ew - octamer_ew.internal_e_ew)


        # ∆G_IGM / kcal mol-1
        dg_igm = Constants.j_mol_to_kcal * (septamer.g + monomer.g - octamer.g)
        dg_igm += dE_opt_to_sp[i-1]

        septamer.calculate_thermochemistry(method='no_s_trans')

        dg_ew = Constants.j_mol_to_kcal * (septamer.g + monomer_ew.g - octamer_ew.g) - TdS_trans + dU_trans
        dg_ew += dE_opt_to_sp[i-1]

        dgs.append(dg_ew)

        print(dg_igm, dg_ew)

    print('∆G   ±(std err)')
    print(np.average(dgs),  np.std(dgs)/(np.sqrt(len(dgs) - 1)))

    # print(f'T∆S = {(s_ig - monomer.s)*temp/(1000 * 4.184)} kcal mol-1')
    # print(f'∆U = {(u_igm - monomer.u) / 1000 / 4.184} kcal mol-1')
