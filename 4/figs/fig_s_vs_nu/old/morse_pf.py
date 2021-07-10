import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from scipy import special
from scipy import integrate
plt.style.use('graph')


def morse_potential(x, d_e, m, nu_hz):
    a = np.sqrt(m * (2.0 * np.pi * nu_hz) ** 2 / (2.0 * d_e))
    return d_e * (1.0 - np.exp(-a * x))**2


def integrand(x, d_e, m, nu_hz, beta):
    morse = morse_potential(x, d_e, m, nu_hz)

    return np.exp(-beta * morse) * special.erf(np.sqrt(beta * (d_e - morse)))


def morse_pf(d_e, nu_hz, kb, t, h, m):

    a = np.sqrt(m * (2.0 * np.pi * nu_hz) ** 2 / (2.0 * d_e))
    x_s = np.linspace((-np.log(2)/a), 5E-9, int(1E4))
    integrands = [integrand(x, d_e, m, nu_hz, (1.0 / (kb * t))) for x in x_s]

    pf = (2.0 / h) * np.sqrt(2.0*  m * np.pi * kb * t) * integrate.trapz(integrands, x_s, 1E5)

    return pf


if __name__ == '__main__':

    kb = 1.38064852E-23     # J K-1
    na = 6.0221409E23       # mol-1
    m_to_ang = 1E10         # A m-1
    c_cm = 2.998E10         # cm s-1
    m = 1.6605E-27          # kg
    d_e = 400 * 1000 / na   # J

    t = 298.15              # K

    nu_cms = np.linspace(100, 2000, 10)
    nu_hzs = nu_cms * c_cm

    min_a = np.sqrt(m * (2.0 * np.pi * min(nu_hzs)) ** 2 / (2.0 * d_e))
    print((-np.log(2)/min_a))
    x_s = np.linspace((-np.log(2)/min_a), 5E-10, int(1E4))
    y_s = np.array([[integrand(x, d_e, m, nu_hz, (1.0 / (kb * t))) for x in x_s] for nu_hz in nu_hzs])
    x_s_ang = m_to_ang * x_s

    norm = mpl.colors.Normalize(vmin=nu_cms.min(), vmax=nu_cms.max())
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
    cmap.set_array([])

    lines = [[(x_s_ang[j], y_s[i][j]) for j in range(len(x_s))] for i in range(len(y_s))]

    line_segments = LineCollection(lines,
                                   linewidths=(0.5, 1, 1.5, 2),
                                   colors=[mpl.cm.jet(x) for x in np.linspace(0.0, 1.0, 10)],
                                   linestyles='solid'
                                   )

    fig, ax = plt.subplots()
    ax.add_collection(line_segments)
    axcb = fig.colorbar(cmap)
    axcb.ax.set_title('$\\nu$ / cm$^{-1}$', fontdict={'fontsize': 10})
    ax.set_xlabel('$x$ / Å')
    ax.set_ylabel('Integrand')
    ax.set_xlim([-2.5, 2.5])
    plt.show()
