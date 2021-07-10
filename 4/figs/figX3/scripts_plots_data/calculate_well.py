import os
import gaptrain as gt
import autode as ade
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
    box_length = np.average(frame.box.size)

    coords = frame.coordinates()

    # Distances from the middle of the box to all the oxygen atoms
    mid_box_dists = np.linalg.norm(coords - box_length / 2.0, axis=1)

    central_oxygen = next(idx for idx in np.argsort(mid_box_dists)
                          if frame.atoms[idx].label == 'O')

    central_o_dists = np.linalg.norm(coords - coords[central_oxygen], axis=1)

    close_o_idxs = [i for i in np.argsort(central_o_dists)
                    if frame.atoms[i].label == 'O'
                    and central_o_dists[i] < max_dist]

    atoms = []
    for o_idx in close_o_idxs:

        # Skip any waters that are on the edge of the box
        if np.max(np.linalg.norm(coords[[o_idx + 1, o_idx + 2]]
                                 - coords[o_idx], axis=1)) > 1.5:
            continue

        # Add the O, and the two hydrogens, which are the next two atoms
        for i in range(3):
            atoms.append(frame.atoms[o_idx + i])

    return ade.Molecule(name='tmp', atoms=atoms)


def get_potential(traj_filename, method, stride=1, max_dist=5):
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

            # Move the central oxygen atom
            for idx in range(3):
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


if __name__ == '__main__':

    n_points = 10

    rs = np.linspace(-2, 2, num=n_points)
    vs = get_potential(traj_filename='dftb_300_nvt.xyz',
                       method=ade.methods.XTB(),
                       stride=8)

    ha_to_kcal = 627.509
    plt.plot(rs, ha_to_kcal * vs, lw=2)
    plt.scatter(rs, ha_to_kcal * vs)

    plt.savefig('potentials')
