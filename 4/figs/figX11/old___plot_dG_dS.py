import matplotlib.pyplot as plt
from scipy import integrate
import numpy as np
plt.style.use('paper')
reds = plt.get_cmap('Reds')
blues = plt.get_cmap('Blues')


class Constants:

    hbar_au = 1.0
    h_au = hbar_au * 2.0 * np.pi
    kb_au = 3.1668114E-6    # hartrees K-1
    h_SI = 6.62607004E-34   # J s

    na = 6.02214086E23       # molecules mol-1
    kb_SI = 1.38064852E-23   # J K-1
    kb_kcalKmol = kb_SI / (4.184 * 1000)  # kcal K-1 mol-1

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
    inverse_ang_inverse_au = 1.0 / 1.88973      # au-1 Å

    k_b = kb_kcalKmol * kcal_mol_to_au          # kcal mol-1 K-1 -> Ha K-1
    r = k_b * n_a                               # Ha K-1


class Molecule:

    @property
    def g_igm(self):
        return (self._g_rrho - self._e) + self._e_sp

    @property
    def g_qrrho_ew(self):
        raise NotImplementedError

    @property
    def g_qrrho(self):
        return (self._g_qrrho - self._e) + self._e_sp

    @property
    def h_igm(self):
        return (self._h_rrho - self._e) + self._e_sp

    @property
    def h_qrrho(self):
        return (self._h_qrrho - self._e) + self._e_sp

    @property
    def de_opt_sp(self):
        return self._e_sp - self._e

    def __init__(self, e, h_qrrho, g_qrrho, e_sp, h_rrho, g_rrho):
        self._e = e
        self._h_qrrho = h_qrrho
        self._g_qrrho = g_qrrho
        self._e_sp = e_sp

        self._h_rrho = h_rrho
        self._g_rrho = g_rrho


class Reaction:

    def delta_energies(self, method='igm', ss='1M'):
        """∆G = x ± y kcal mol-1"""
        dgs = np.array([627.509*(getattr(self.items[i], f'g_{method}')
                                 + getattr(self.items[0], f'g_{method}')
                                 - getattr(self.items[1], f'g_{method}'))
                        for i in range(2, 10)])

        dhs = np.array([627.509*(getattr(self.items[i], f'h_{method}')
                                 + getattr(self.items[0], f'h_{method}')
                                 - getattr(self.items[1], f'h_{method}'))
                        for i in range(2, 10)])

        dsts = dgs - dhs   # -T∆S / kcal mol-1

        rt = (8.3145 * 298.15) / (1000 * 4.184)    # kcal mol-1

        if ss == '1M':
            correction_kcal = rt * np.log(40605/1661)

        elif ss == 'liq':
            def eff_vol(molarity):
                """Effective volume in molecules Å-3"""
                return molarity * 0.001 * 1E30 / 6.0221409E23

            correction_kcal = rt * (np.log(40605 / eff_vol(molarity=55.5/7))
                                    + np.log(40605 / eff_vol(molarity=55.5))
                                    - np.log(40605 / eff_vol(molarity=55.5/8)))

        else:
            raise NotImplementedError

        dgs += correction_kcal
        dsts += correction_kcal

        return (np.mean(dgs),
                np.std(dgs)/(np.sqrt(len(dgs) - 1)),

                np.mean(dsts),
                np.std(dsts)/(np.sqrt(len(dsts) - 1)))

    def __init__(self, *args):
        """
        Reaction in order

        monomer
        octamer
        septamer1
        septamer2
        septamer3
        septamer4
        septamer5
        septamer6
        septamer7
        septamer8

        :param args:
        """
        self.items = args


def plot_dG():
    ax.bar(x=[0, 0.5, 1.5, 2.0, 3.0, 3.5],
           width=0.5,
           height=[e[0] for e in delta_es],
           yerr=[e[1] for e in delta_es],
           linewidth=0.5,
           color=[blues(0.4), reds(0.4),
                  blues(0.5), reds(0.5),
                  blues(0.6), reds(0.6)],
           edgecolor='darkgray',
           capsize=2)

    ax.set_ylabel('∆$G$ / kcal mol$^{-1}$')


def plot_dS():
    ax.bar(x=[0, 0.5, 1.5, 2.0, 3.0, 3.5],
           width=0.5,
           height=[e[2] for e in delta_es],
           yerr=[e[3] for e in delta_es],
           linewidth=0.5,
           color=[blues(0.4), reds(0.4),
                  blues(0.5), reds(0.5),
                  blues(0.6), reds(0.6)],
           hatch='/',
           capsize=2)

    ax.yaxis.set_ticks_position("right")
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('-T∆$S$ / kcal mol$^{-1}$')
    return


if __name__ == '__main__':

    smd = Reaction(Molecule(-76.288703443144, -76.26360045460233, -76.28514535940438,  -76.350985605823 , -76.26360045460233, -76.28514535919648),
                   Molecule(-610.433983658315, -610.2140982219828, -610.277280121547,  -610.840350617938, -610.2140982219828, -610.281451708829),
                   Molecule(-534.12940950742, -533.9371068609033, -533.9942137795163,  -534.483816573074, -533.9371068609033, -533.996445983519),
                   Molecule(-534.12912670771, -533.9377740081061, -533.9935062669732,  -534.484124253745, -533.9377740081061, -533.9953996468886),
                   Molecule(-534.126645015954, -533.935632927275, -533.9921377225431,  -534.484446846104, -533.935632927275, -533.9948616348487),
                   Molecule(-534.123952514425, -533.9321026652082, -533.9909755630321, -534.483679411155, -533.9321026652082, -533.9952296133653),
                   Molecule(-534.122528497788, -533.9317335339696, -533.9900214595227, -534.482666928205, -533.9317335339696, -533.9941873317987),
                   Molecule(-534.129192138055, -533.9381699483047, -533.9939597762663, -534.485678939561, -533.9381699483047, -533.9964300553619),
                   Molecule(-534.121709571807, -533.9317638621948, -533.9872949274167, -534.479553024392, -533.9317638621948, -533.9893161318665),
                   Molecule(-534.128584701367, -533.9374692734386, -533.9935968301443, -534.484199274648, -533.9374692734386, -533.9962259260279))

    for i in range(8):
        print('∆E_opt->sp =',
              smd.items[i+2].de_opt_sp + smd.items[0].de_opt_sp - smd.items[1].de_opt_sp)

    exit()
    # ∆G ±
    print(smd.delta_energies(method='qrrho', ss='liq')[:2])
    exit()

    # -------------------------- Plot --------------------------------
    fig, ax = plt.subplots()

    delta_es = [smd.delta_energies(method='igm', ss='1atm'),
                smd.delta_energies(method='igm', ss='liq'),
                smd.delta_energies(method='qrrho', ss='liq')]

    # plot_dG()
    # plot_dS()

    ax.plot([-0.5, 4], [0, 0], c='k', lw=1.5)

    ax.set_xticks([])
    plt.xlim(-0.5, 4)
    plt.ylim(-12, 10)

    plt.tight_layout()
    plt.savefig('water_in_water_TdS.pdf')
