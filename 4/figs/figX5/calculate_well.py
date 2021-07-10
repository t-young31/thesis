import os
import gaptrain as gt
import autode as ade
from autode.atoms import Atom
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('paper')
ade.Config.n_cores = 4


def get_truncated_portion(frame, max_dist):
    """
    Given a frame discard all but the central water molecule and the water
    molecules that are within max_dist of it

    :param frame:
    :param max_dist:
    :return:
    """
    a, b, c, = frame.box.size
    solvent_n_atoms = system.molecules[-1].n_atoms

    atom = frame.atoms.pop(133)
    frame.atoms.insert(49, atom)

    atom = frame.atoms.pop(134)
    frame.atoms.insert(50, atom)

    atom = frame.atoms.pop(134)
    frame.atoms.insert(118, atom)

    n_atoms = len(frame.atoms)

    # Build a supercell 3x3x3
    ss_coords = [frame.coordinates()]   # Tensor 27 x n_atoms x 3

    for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for k in [-1, 0, 1]:

                if i == j == k == 0:
                    continue

                coords = frame.coordinates()
                coords[:, 0] += i * a
                coords[:, 1] += j * b
                coords[:, 2] += k * c

                ss_coords.append(coords)

    ss_coords = np.array(ss_coords)

    s_coord = ss_coords[0, 0, :]     # Position of the first atom in the solute

    # Select the nearest solute atoms to the first one, in the central box
    # of the supercell
    atoms = []
    for idx in range(solute.n_atoms):
        dists = np.linalg.norm(ss_coords[:, idx] - s_coord, axis=1)
        coord = ss_coords[np.argmin(dists), idx, :]
        atoms.append(Atom(frame.atoms[idx].label, *coord))

    # Now grab the first atoms for every solvent molecule within max_dist
    flat_ss_coords = ss_coords.reshape(27*n_atoms, 3)

    # This took agesss to work out the indexing for :(
    dists = np.linalg.norm(flat_ss_coords - s_coord, axis=1).reshape(27, n_atoms)
    min_idxs = np.argmin(dists, axis=0)

    for i in range(solute.n_atoms, n_atoms)[::solvent_n_atoms]:

        # If the first atom of this solute molecule is outside the cut-off
        # radius
        cell = min_idxs[i]
        if dists[cell, i] > max_dist:
            continue

        # Coordinate of the first atom in the solvent
        s_coord = ss_coords[cell, i, :]

        # now add all the solvent atoms that are closest to this one
        # out of the full supercell
        for j in range(solvent_n_atoms):
            s_dists = np.linalg.norm(ss_coords[:, i+j] - s_coord, axis=1)
            coord = ss_coords[np.argmin(s_dists), i+j, :]

            atoms.append(Atom(frame.atoms[i+j].label, *coord))

    return ade.Molecule(name='tmp', atoms=atoms)


def get_potential(traj_filename, method, stride=1, max_dist=7):
    """

    :param traj_filename:
    :param method:
    :param stride:
    :param max_dist:
    :return:
    """
    data_filename = f'{traj_filename}_potentials.txt'

    if os.path.exists(data_filename):
        all_potentials = np.loadtxt(data_filename)
        return np.average(all_potentials, axis=0)

    frames = gt.Data(traj_filename)
    all_potentials = []   # arrays of V(r) for each frame in the trajectory

    for i, frame in enumerate(frames[::stride]):

        mol = get_truncated_portion(frame, max_dist=max_dist)
        potentials = []

        for r in rs:
            shift_mol = mol.copy()
            shift_mol.name = f'{i}_{r:.3f}'

            # Move the the solute
            for idx in range(solute.n_atoms):
                shift_mol.atoms[idx].translate(vec=np.array([r, 0.0, 0.0]))

            shift_mol.single_point(method=method)

            potentials.append(shift_mol.energy)

            # Clean up all the generated files
            for filename in os.listdir(os.getcwd()):
                if filename.startswith(shift_mol.name):
                    os.remove(filename)

        # Append the relative energies
        potentials = np.array(potentials) - min(potentials)
        all_potentials.append(potentials)

    all_potentials = np.array(all_potentials)
    np.savetxt(data_filename, all_potentials)

    return np.average(all_potentials, axis=0)


def rough_plot():
    """Generate a rough plot of V vs x for quick visualisation"""

    ha_to_kcal = 627.509
    plt.plot(rs, ha_to_kcal * vs, lw=2)
    plt.scatter(rs, ha_to_kcal * vs)

    plt.savefig('potentials')
    return None


def generate_trajectory():
    """Generate a short DFTB trajectory"""

    config = system.random(min_dist_threshold=1.6)
    config.save()
    exit()
    config.run_dftb(max_force=0.1)

    # Short 10 ps equilibration
    eqm = gt.md.run_dftbmd(config,
                           temp=300,
                           dt=0.5,
                           interval=10,
                           ps=10)

    # And a slightly longer NVT simulation
    nvt = gt.md.run_dftbmd(eqm[-1],
                           temp=300,
                           dt=0.5,
                           interval=10,
                           ps=20)
    for config in nvt:
        config.wrap()

    nvt.save(filename='dftb_300_nvt.xyz')
    return None


if __name__ == '__main__':

    solute = gt.Molecule('CO2.xyz')
    solvent_mw = 41.05

    system = gt.System(solute, box_size=[12, 12, 12])

    # Volume available for the solvent, in Å^3
    solvent_v = np.product(system.box.size) - (4*np.pi/3 * solute.radius**3)
    # and the density, in molecules Å-3
    solvent_p = ((0.786 * 1E6 / solvent_mw) * 6.022E23 / 1E30)

    # Add the correct number of solvent molecules to get to the correct solvent
    # density
    system.add_molecules(gt.Molecule('C6H6.xyz'),
                         n=int(solvent_v * solvent_p))

    # generate_trajectory()

    n_points = 13
    rs = np.linspace(-2.5, 2.5, num=n_points)
    vs = get_potential(traj_filename='dftb_300_nvt.xyz',
                       method=ade.methods.XTB(),
                       stride=8)
    rough_plot()
