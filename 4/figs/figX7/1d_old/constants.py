import numpy as np


class Constants(object):

    hbar_au = 1.0
    h_au = hbar_au * 2.0 * np.pi
    kb_au = 3.1668114E-6    # hartrees K-1
    h_SI = 6.62607004E-34   # J s

    na = 6.02214086E23       # molecules mol-1
    kb_SI = 1.38064852E-23   # J K-1
    kb_JKmol = kb_SI * na    # J K-1 mol-1

    k_b = 1.38064852E-23                        # J K-1
    h = 6.62607004E-34                          # J s
    n_a = 6.022140857E23                        # molecules mol-1
    atm_to_pa = 101325                          # Pa
    dm_to_m = 0.1                               # m
    amu_to_kg = 1.660539040E-27                 # Kg
    r = k_b * n_a                               # J K-1 mol-1
    c = 299792458                               # m s-1
    c_in_cm = c * 100                           # cm s-1
    ang_to_m = 1E-10                            # m
    m_to_ang = 1E10                             # Å
    m_to_bohr = 1.89e+10                        # au m-1
    amu_to_au = 1822.888486                     # m_e amu-1
    kj_mol_to_au = 0.00038087980                # Ha (kJ mol-1)-1       Hartrees = au of energy
    inverse_ang_inverse_au = 1.0 / 1.88973      # au-1 Å
