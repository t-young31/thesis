import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import special
plt.style.use('paper')

c_cm = 2.998E10         # cm s-1
m_to_ang = 1E10         # A m-1

kb = 1.38064852E-23     # J K-1
na = 6.0221409E23       # mol-1
r = 8.3144598 / (1000 * 4.184)          # kcal K-1 mol-1
h = 6.62607004E-34      # J s
T = 298.15              # K
beta = 1.0 / (kb * T)   # J-1

r_e = 1E-10             # 1 Å
m = 1.6605E-27          # kg

d_e = 400 * 1000 / na   # J


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

    pf = ((2.0 / h)
          * np.sqrt(2.0 * m * np.pi * kb * t)
          * np.trapz(integrands, x_s))

    return pf


def S_HO_classical(nu_per_cm):
    """Classical entropy for a 1D ossilator with ν = nu_per_cm cm^-1"""
    return r * (np.log(np.exp(1) / (beta * h * nu_per_cm * c_cm)))


def S_HO_qunatum(nu_per_cm):
    """Quantum entropy for a 1D ossilator with ν = nu_per_cm cm^-1"""

    h_nu_beta = beta * h * (nu_per_cm * c_cm)
    s = (r * (h_nu_beta / (np.exp(h_nu_beta) - 1.0)
              - np.log(1.0 - np.exp(- h_nu_beta))))

    return s


def S_truncatedHO_classical(nu_per_cm):

    def z(nu):
        return ((np.pi * kb * T / (4.0 * h * nu)) * special.erfc(-r_e * np.sqrt(beta / (2.0 * m))) *
                special.erfc(-r_e * np.sqrt(beta * m * (2.0 * np.pi * nu)**2 / (2.0 * m))))          # nu = omega / 2 pi

    def dZ_dT(nu):
        a = r_e**2 / (2.0 * kb * m)
        b = r_e**2 * m * (2.0 * np.pi * nu)**2 / (2.0 * kb)
        term0 = special.erfc(-np.sqrt(a/T)) * special.erfc(-np.sqrt(b/T))
        term1 = (b * np.exp(-b / T) * special.erfc(-np.sqrt(a/T))) / (np.sqrt(b * np.pi * T))
        term2 = (a * np.exp(-a / T) * special.erfc(-np.sqrt(b/T))) / (np.sqrt(a * np.pi * T))

        return (kb / (4.0 * h * nu)) * (term0 + term1 + term2)

    pf = z(nu=nu_per_cm * c_cm)
    return r * (T * (1.0 / pf) * dZ_dT(nu_per_cm * c_cm) + np.log(pf))


def S_morse_classical(nu_per_cm):

    def dZ_dT(nu):

        dt = 1E-8
        return (morse_pf(d_e, nu, kb, T + dt, h, m)
                - morse_pf(d_e, nu, kb, T, h, m)) / dt

    pf = morse_pf(d_e, nu_per_cm*c_cm, kb, T, h, m)

    return r * (T * (1.0 / pf) * dZ_dT(nu_per_cm*c_cm) + np.log(pf))


def plot_HO_quant_class():

    fig, ax1 = plt.subplots()

    ax1.plot(nus, T * S_HO_qunatum(nus), label='quantum', lw=1.5)
    ax1.plot(nus, T * S_HO_classical(nus), label='classical', lw=1.5)
    ax1.plot([-100, 3000], [0.0, 0.0], c='k', ls='--')
    ax1.set_xlim(-20, 3000)
    ax1.set_ylim(-1.2, 5)
    ax1.set_xlabel('$\\nu$ / cm$^{-1}$')
    ax1.set_ylabel('$TS_{HO}$ / kcal mol$^{-1}$')
    ax1.legend()

    plt.tight_layout()
    plt.savefig('s_vs_n.pdf')

    return None


def plot_classicals():

    fig, ax2 = plt.subplots()

    ax2.plot(nus, T * S_HO_classical(nus),
             label='HO', color='tab:orange', lw=1.5)

    ax2.plot(nus, T * S_truncatedHO_classical(nus),
             label='trunc. HO', color='tab:green', lw=1.5)

    if not os.path.exists('morse_TS.txt'):
        morse_vals = np.array([T * S_morse_classical(nu) for nu in nus])
        np.savetxt('morse_TS.txt', morse_vals)

    ax2.plot(nus, np.loadtxt('morse_TS.txt'),
             label='Morse', color='tab:red', lw=1.5)

    ax2.plot([-100, 3000], [0.0, 0.0], c='k', ls='--')
    ax2.set_xlim(-20, 3000)
    ax2.set_ylim(-1.2, 5)
    ax2.set_xlabel('$\\nu$ / cm$^{-1}$')
    ax2.set_ylabel('$TS_{classical}$ / kcal mol$^{-1}$')
    ax2.legend()

    plt.tight_layout()
    plt.savefig('s_vs_nu_trunc_morse.pdf')

    return None


if __name__ == '__main__':

    nus = np.linspace(0.1, 3000, num=500)

    plot_classicals()
    plot_HO_quant_class()
