import autode as ade
import numpy as np
from autode.methods import ORCA
from autode.mol_graphs import make_graph
from autode.pes.pes_2d import PES2d
ade.Config.n_cores = 24
n_points = 10

if __name__ == '__main__':

    reac = ade.Reactant('sn2_init.xyz', charge=-1)
    make_graph(reac)

    prod = ade.Product('sn2_final.xyz', charge=-1)
    make_graph(prod)

    pes = PES2d(reac, prod,
                r1s=np.linspace(3.4, 1.3, n_points),
                r1_idxs=(0, 2),                          # F-C
                r2s=np.linspace(1.7, 2.9, n_points),
                r2_idxs=(2, 1))                          # C-Cl
    pes.calculate(name='orca_sn2_surface',
                  method=ORCA(),
                  keywords=ade.OptKeywords(['PBE0', 'ma-def2-SVP', 'LooseOpt']))

    energies = np.zeros(shape=(n_points, n_points))
    for i in range(n_points):
        for j in range(n_points):
            energies[i, j] = pes.species[i, j].energy

    np.savetxt('orca_sn2_surface.txt', energies)
