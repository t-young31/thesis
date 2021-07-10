import autode as ade
from autode.neb import NEB
ade.Config.n_cores = 8
ade.Config.adaptive_neb_k = True

 
if __name__ == '__main__':

    r, p = ade.Reactant('r.xyz', mult=2), ade.Product('p.xyz', mult=2)
    for mol in (r, p):
        mol.optimise(method=ade.methods.ORCA())
        mol.print_xyz_file()

    neb = NEB(initial_species=r,
              final_species=p,
              num=6)

    neb.calculate(method=ade.methods.ORCA(),
                  n_cores=10)
