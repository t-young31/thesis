import numpy as np
import autode as ade
from autode.atoms import get_vdw_radius


if __name__ == '__main__':

    h2o = ade.Molecule(smiles='O')
    coords = h2o.coordinates
    centroid = np.average(coords, axis=0)
    r_o = np.linalg.norm(coords[0] - centroid) + get_vdw_radius('O')
    r_h = np.linalg.norm(coords[0] - centroid) + get_vdw_radius('H')

    r = max(r_h, r_o)
    v_m = (4.0 * np.pi / 3.0) * r**3
    print(v_m)

    v_free = ((18.02 / 1000)         # M_w     in kg mol-1
              /
              (6.022E23             # N_a     in mol-1
               * 1000)              # density in kg m^3
              ) * 1E30              # m^3 -> Å^3

    v_free -= v_m

    N_c = 1
    v_c = (v_m**(1/3) + v_free**(1/3))**3

    V = N_c * v_c
    print(V)
    print(f'l = {V**(1/3)}')
