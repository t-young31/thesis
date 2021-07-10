import os
import sys
import gaptrain as gt
import autode as ade
import numpy as np
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


def save_potential(traj_filename, stride=1, max_dist=5, frame_n=0):
    """

    :param traj_filename:
    :param stride:
    :param max_dist:
    :return:
    """
    frames = gt.Data(traj_filename)

    for i, frame in enumerate(frames[::stride]):

        if i != frame_n:
            continue

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

        # Save the potential array to a .txt
        potentials = np.array(potentials) - min(potentials)
        np.savetxt(f'{traj_filename}_{frame_n}_pbe0_potential.txt',
                   potentials)

    return


if __name__ == '__main__':

    method = ade.methods.ORCA()
    method.keywords.sp.basis_set = 'def2-SVP'

    n_points = 7

    rs = np.linspace(-1.5, 1.5, num=n_points)
    save_potential(traj_filename='GAP_rPBE0-D3_nvt_300K.xyz',
                   stride=8,
                   frame_n=int(sys.argv[1]))

