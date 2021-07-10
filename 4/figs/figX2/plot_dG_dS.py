import matplotlib.pyplot as plt
import numpy as np
plt.style.use('paper')
reds = plt.get_cmap('Reds')
blues = plt.get_cmap('Blues')


class Molecule:

    @property
    def g_igm(self):
        return (self._g_rrho - self._e) + self._e_sp

    @property
    def g_qrrho(self):
        return (self._g_qrrho - self._e) + self._e_sp

    @property
    def h_igm(self):
        return (self._h_rrho - self._e) + self._e_sp

    @property
    def h_qrrho(self):
        return (self._h_qrrho - self._e) + self._e_sp

    def __init__(self, e, h_qrrho, g_qrrho, e_sp, h_rrho, g_rrho):
        self._e = e
        self._h_qrrho = h_qrrho
        self._g_qrrho = g_qrrho
        self._e_sp = e_sp

        self._h_rrho = h_rrho
        self._g_rrho = g_rrho


class Reaction:

    def delta_energies(self, method='igm', ss='1atm'):
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

        if ss == '1atm':
            # Entropies already calculated at 1 atm
            correction_kcal = 0

        # elif ss == '1M':
        #    correction_kcal = rt * np.log(40605/1661)

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

        return np.mean(dgs), np.std(dgs)/(np.sqrt(len(dgs) - 1)), np.mean(dsts), np.std(dsts)/(np.sqrt(len(dsts) - 1))

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


if __name__ == '__main__':

    # Gas phase PBE0-D3BJ/def2-SVP loosopt // DLPNO-CCSD(T)/def2-TZVPP at 1atm
    gas = Reaction(Molecule(-76.276608932034,-76.25124030960154,-76.27277719424667, -76.336829911324  , -76.25124030960154, -76.27277719419038),
                   Molecule( -610.394000579583,-610.170453206054,-610.2313972386472,-610.793529928081 , -610.170453206054, -610.2347361697064),
                   Molecule(-534.09564374938,-533.8999961068107,-533.9551116640462,-534.444328310093  , -533.8999961068107, -533.9568583888058),
                   Molecule(-534.092518621825,-533.8968908613684,-533.9523145899593,-534.441451566438 , -533.8968908613684, -533.954189425665),
                   Molecule( -534.093552256922,-533.8973363078736,-533.9529510348643,-534.441562209708, -533.8973363078736, -533.9545842149414),
                   Molecule( -534.087280639382,-533.8932220797443,-533.9483432033645,-534.440658569016, -533.8932220797443, -533.9508978010879),
                   Molecule( -534.086946981991,-533.8926382107583,-533.9478045884676,-534.440014006283, -533.8926382107583, -533.9499956462853),
                   Molecule( -534.094461884023,-533.8990913065269,-533.9549580384031,-534.444607770227, -533.8990913065269, -533.957552050183),
                   Molecule( -534.096093703941,-533.9000220604281,-533.9548599134145,-534.444475376512, -533.9000220604281, -533.9566015289574),
                   Molecule( -534.09235340845,-533.8967925745227,-533.9523585588456,-534.441371188255 , -533.8967925745227, -533.954219144212))

    print('(∆E_TZ - ∆E_QZ)=',
          627.5*(-76.361116805936 -534.610399361119 --610.983263498505) -
          627.5*(-76.336829911324 -534.444475376512 --610.793529928081))

    ss = 627.509*(-76.28212854985226 - -76.28514535940438)
    print(f'1 atm -> 1 M = {ss:.5f}')

    cpcm = Reaction(
        Molecule(-76.287389105275, -76.2622307402399, -76.28377435710208,   -76.347693642677 , -76.2622307402399, -76.28377435692312),
        Molecule(-610.43101614605, -610.209946654115, -610.2717100917978,   -610.829285184708, -610.209946654115, -610.2750579386827),
        Molecule(-534.122598976861, -533.9295606840542, -533.9853769410822, -534.472642253303, -533.9295606840542, -533.9875709328695),
        Molecule(-534.117343503336, -533.9239691528271, -533.9825088171939, -534.469554078253, -533.9239691528271, -533.9852982354548),
        Molecule(-534.124014622961, -533.9321592967783, -533.9875826589193, -534.474750119401, -533.9321592967783, -533.9897978744827),
        Molecule(-534.121057830409, -533.9282868245107, -533.9860036755525, -534.473365830876, -533.9282868245107, -533.9894405750998),
        Molecule(-534.119258624784, -533.9274209224809, -533.9836697924085, -534.472193144436, -533.9274209224809, -533.9862661961873),
        Molecule(-534.12650346481, -533.9336197722118, -533.990247966017,   -534.475813606831, -533.9336197722118, -533.9929946691295),
        Molecule(-534.120911166469, -533.9278727742018, -533.985314325824,  -534.470640350527, -533.9278727742018, -533.9879248711704),
        Molecule(-534.125168899498, -533.9330243196298, -533.9879429966388, -534.474373198342, -533.9330243196298, -533.9899660594167)
    )

    smd = Reaction(
        Molecule(-76.288703443144, -76.26360045460233, -76.28514535940438,  -76.350985605823 , -76.26360045460233, -76.28514535919648),
        Molecule(-610.433983658315, -610.2140982219828, -610.277280121547,  -610.840350617938, -610.2140982219828, -610.281451708829),
        Molecule(-534.12940950742, -533.9371068609033, -533.9942137795163,  -534.483816573074, -533.9371068609033, -533.996445983519),
        Molecule(-534.12912670771, -533.9377740081061, -533.9935062669732,  -534.484124253745, -533.9377740081061, -533.9953996468886),
        Molecule(-534.126645015954, -533.935632927275, -533.9921377225431,  -534.484446846104, -533.935632927275, -533.9948616348487),
        Molecule(-534.123952514425, -533.9321026652082, -533.9909755630321, -534.483679411155, -533.9321026652082, -533.9952296133653),
        Molecule(-534.122528497788, -533.9317335339696, -533.9900214595227, -534.482666928205, -533.9317335339696, -533.9941873317987),
        Molecule(-534.129192138055, -533.9381699483047, -533.9939597762663, -534.485678939561, -533.9381699483047, -533.9964300553619),
        Molecule(-534.121709571807, -533.9317638621948, -533.9872949274167, -534.479553024392, -533.9317638621948, -533.9893161318665),
        Molecule(-534.128584701367, -533.9374692734386, -533.9935968301443, -534.484199274648, -533.9374692734386, -533.9962259260279)
    )

    # -------------------------- Plot --------------------------------
    fig, ax = plt.subplots()

    delta_es = [cpcm.delta_energies(method='igm', ss='1atm'),
                smd.delta_energies(method='igm', ss='1atm'),

                cpcm.delta_energies(method='igm', ss='liq'),
                smd.delta_energies(method='igm', ss='liq'),

                cpcm.delta_energies(method='qrrho', ss='liq'),
                smd.delta_energies(method='qrrho', ss='liq')
                ]

    if False:
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

    if True:
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

    ax.plot([-0.5, 4], [0, 0], c='k', lw=1.5)

    ax.set_xticks([])
    plt.xlim(-0.5, 4)
    plt.ylim(-12, 10)

    plt.tight_layout()
    plt.savefig('water_in_water_TdS.pdf')
