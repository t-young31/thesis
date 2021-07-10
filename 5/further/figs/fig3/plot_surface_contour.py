from scipy.interpolate import interp2d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import autode as ade
from autode.constants import Constants
plt.style.use('paper')
n_points = 10


def drs_to_abs(idxs, drs):
    r = init_point.distance(*idxs)
    rs = []

    for dr in drs:
        r += dr
        rs.append(r)

    return rs


if __name__ == '__main__':

    fig, ax = plt.subplots()

    energies = np.loadtxt('orca_sn2_surface.txt')
    energies -= np.min(energies)
    energies *= Constants.ha2kcalmol

    x = np.linspace(3.4, 1.3, n_points)
    y = np.linspace(1.7, 2.9, n_points)

    init_point = ade.Molecule('path_opt0_orca.xyz')

    drs_CF_min0p05 = [-0.2999963166993746, -0.299947664236493, -0.3, -0.3, -0.3, -0.3, 0.0]
    drs_CCl_min0p05 = [0.08428161307157019, 0.07045741945647686, 0.07885138210871159, 0.1315327607473873, 0.2857430136680222, 0.3, 0.3]

    ax.plot(drs_to_abs(idxs=(1, 2), drs=drs_CCl_min0p05),
            drs_to_abs(idxs=(0, 2), drs=drs_CF_min0p05),
            marker='o', c='b',
            label='adapt. $\Delta r_{min} = 0.05$')

    drs_CCl_min0p01 = [0.0499407206740569, 0.041561634312418556, 0.053505246720705814, 0.11571772333026374, 0.29112079728208073, 0.3, 0.3]
    drs_CF_min0p01 = [-0.2999991997626649, -0.2996393063089389, -0.3, -0.29934000393180327, -0.3, -0.3, 0.0]

    ax.plot(drs_to_abs(idxs=(1, 2), drs=drs_CCl_min0p01),
            drs_to_abs(idxs=(0, 2), drs=drs_CF_min0p01),
            marker='o', c='m',
            label='adapt. $\Delta r_{min} = 0.01$')

    # Linear path
    ax.plot([init_point.distance(1, 2) + i * (2.85 - init_point.distance(1, 2))/10 for i in range(11)],
            [init_point.distance(0, 2) -0.05 + i * (1.4 - init_point.distance(0, 2))/10 for i in range(11)],
            marker='o', c='r', label='linear')

    rs_CF_irc = [2.0578, 2.1523, 2.255, 2.3449, 2.4398, 2.5264, 2.6107, 2.6848, 2.6962, 2.7078, 2.7257, 2.745, 2.765, 2.7838, 2.804, 2.8216, 2.8411, 2.8594, 2.8776, 2.8957][::-1] + [1.9277, 1.8295, 1.6695, 1.5781, 1.4767, 1.4009, 1.4021, 1.3991, 1.3959, 1.4035, 1.3996]
    rs_CCl_irc = [2.1394, 2.0443, 1.9504, 1.8747, 1.8191, 1.8157, 1.8085, 1.8106, 1.8088, 1.8078, 1.8068, 1.8065, 1.8055, 1.8055, 1.8039, 1.8049, 1.8032, 1.803, 1.8028, 1.8026][::-1] + [2.2548, 2.3483, 2.4889, 2.5699, 2.6602, 2.7485, 2.8275, 2.9144, 2.9925, 3.0713, 3.0806]

    ax.plot(rs_CCl_irc, rs_CF_irc,
            c='green', label='IRC', marker='s', ms=3)

    ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%0.1f'))

    plt.scatter([1.993], [2.197], c='k', s=40, label='TS', zorder=20)

    im = plt.imshow(energies,
                    extent=[1.7, 2.9, 1.3, 3.4],
                    cmap='RdBu_r',
                    alpha=0.6,
                    aspect=0.55,
                    vmax=100
                    )

    cbar = plt.colorbar(im)
    cbar.set_label('∆$E$ / kcal mol$^{-1}$')
    # plt.xlim(1.5, 3.1)
    # plt.ylim(1.5, 3.1)
    plt.ylabel('$r_{C-F}$ / Å')
    plt.xlabel('$r_{C-Cl}$ / Å')
    plt.legend(framealpha=0.5, fontsize=9)
    plt.tight_layout()
    plt.savefig('adapt_surface_sn2.png', dpi=400)
