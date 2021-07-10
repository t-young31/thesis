import gaptrain as gt
gt.GTConfig.n_cores = 8


if __name__ == '__main__':

    system = gt.System(box_size=[10, 10, 10])
    system.add_solvent('h2o', n=33)

    config = system.random()
    config.run_dftb(max_force=0.1)

    eqm = gt.md.run_dftbmd(config,
                           temp=300,
                           dt=0.5,
                           interval=10,
                           ps=10)

    nvt = gt.md.run_dftbmd(eqm[-1],
                           temp=300,
                           dt=0.5,
                           interval=10,
                           ps=20)
    for config in nvt:
        config.wrap()

    nvt.save(filename='dftb_300_nvt.xyz')
