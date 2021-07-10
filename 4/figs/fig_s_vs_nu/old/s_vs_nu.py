from morse_pf import *
plt.style.use('paper')


def morse_potential(r, d_e, nu, r_e ):
    a = np.sqrt(m * (2.0 * np.pi * nu) ** 2 / (2.0 * d_e))
    return d_e * (1.0 - np.exp(-a * (r - r_e)))**2


def ho_potential(r, nu, r_e):
    return 0.5 * m * (2.0 * np.pi * nu)**2 * (r - r_e)**2


def qunatum(nus_hz):

    entropies = []

    for nu in nus_hz:
        h_nu_over_kb_t = (h * nu) / (kb * t)
        s = r * (h_nu_over_kb_t / (np.exp(h_nu_over_kb_t) - 1.0) - np.log(1.0 - np.exp(- h_nu_over_kb_t)))
        entropies.append(s)

    return plt.plot(nus, entropies, label='quantum', c='#1f77b4')


def classical(nus_hz):

    classical_entropies = []

    for nu in nus_hz:
        h_nu_over_kb_t = (h * nu) / (kb * t)
        s_clasical = r * (np.log(np.exp(1) / h_nu_over_kb_t))
        classical_entropies.append(s_clasical)

    return plt.plot(nus, classical_entropies, label='HO', c='#ff7f0e')


def truncated_classical(nus_hz):

    classical_entropies = []

    def z(nu):
        return ((np.pi * kb * t / (4.0 * h * nu)) * special.erfc(-r_e * np.sqrt(beta / (2.0 * m))) *
                special.erfc(-r_e * np.sqrt(beta * m * (2.0 * np.pi * nu)**2 / (2.0 * m))))          # nu = omega / 2 pi

    def dZ_dT(nu):
        a = r_e**2 / (2.0 * kb * m)
        b = r_e**2 * m * (2.0 * np.pi * nu)**2 / (2.0 * kb)
        term0 = special.erfc(-np.sqrt(a/t)) * special.erfc(-np.sqrt(b/t))
        term1 = (b * np.exp(-b / t) * special.erfc(-np.sqrt(a/t))) / (np.sqrt(b * np.pi * t))
        term2 = (a * np.exp(-a / t) * special.erfc(-np.sqrt(b/t))) / (np.sqrt(a * np.pi * t))

        return (kb / (4.0 * h * nu)) * (term0 + term1 + term2)

    for nu in nus_hz:
        pf = z(nu)
        s_clasical = r * (t * (1.0 / pf) * dZ_dT(nu) + np.log(pf))
        classical_entropies.append(s_clasical)

    return plt.plot(nus, classical_entropies, label='trunc. HO', c='#2ca02c')


def morse_classical(nus_hz):

    classical_entropies = []

    def dZ_dT(nu, temp):

        dt = 1E-8
        return (morse_pf(d_e, nu, kb, temp + dt, h, m) - morse_pf(d_e, nu, kb, temp, h, m)) / dt

    for nu in nus_hz:
        pf = morse_pf(d_e, nu, kb, t, h, m)
        s_clasical = r * (t * (1.0 / pf) * dZ_dT(nu, t) + np.log(pf))
        classical_entropies.append(s_clasical)

    return plt.plot(nus, classical_entropies, label='Morse', c='#d62728')


def plot(nus_hz):
    # qunatum(nus_hz)
    classical(nus_hz)
    truncated_classical(nus_hz)
    morse_classical(nus_hz)

    plt.plot(np.linspace(-200, 5000, 10), np.zeros(10), ls='--', c='k')
    plt.xlabel('$\\nu$ / cm$^{-1}$')
    plt.ylabel('$S$ / J K$^{-1}$ mol$^{-1}$')
    plt.ylim(-20)
    plt.xlim([-22, 3000])
    plt.legend()

    return plt.show()


if __name__ == '__main__':

    nus = np.linspace(0.1, 3000, 100)
    c_cm = 2.998E10         # cm s-1
    m_to_ang = 1E10         # A m-1
    nus_hz = nus * c_cm

    kb = 1.38064852E-23     # J K-1
    na = 6.0221409E23       # mol-1
    r = 8.3144598           # J K-1 mol-1
    h = 6.62607004E-34      # J s
    t = 298.15              # K
    beta = 1.0 / (kb * t)   # J-1

    r_e = 1E-10             # 1 Å
    m = 1.6605E-27          # kg

    d_e = 400 * 1000 / na   # J

    # plot(nus_hz)

    nu = 2000 * c_cm
    xs = np.linspace(-1.5E-10, 4E-10, 500)
    ys_ho = np.array([ho_potential(x, nu, r_e=0.0) for x in xs])
    ys_ho_trunc = np.array([ho_potential(x, nu, r_e=0.0) if x > -r_e else 1E10 for x in xs])
    ys_morse = np.array([morse_potential(x, d_e, nu, r_e=0.0) for x in xs])

    plt.plot(xs * m_to_ang, 1E18 * ys_ho,  c='#ff7f0e', ls='--')
    plt.plot(xs * m_to_ang, 1E18 * ys_ho_trunc,  c='#2ca02c')
    plt.xlabel('x / Å')
    plt.ylabel('V(x) / $\\times 10^{-18}\;$ J')
    plt.ylim([0.0, 2.0])
    plt.show()
    plt.plot(xs * m_to_ang, 1E18 * ys_ho,  c='#ff7f0e', ls='--')
    plt.plot(xs * m_to_ang, 1E18 * ys_morse, c='#d62728')
    plt.xlabel('x / Å')
    plt.ylabel('V(x) / $\\times 10^{-18}\;$ J')
    plt.ylim([0.0, 2.0])
    plt.show()
