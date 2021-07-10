import gaptrain as gt
import numpy as np
import matplotlib.pyplot as plt
blues = plt.get_cmap('Blues')
plt.style.use('paper')


def plot_line(energies, traj, label, color):

    dists = []
    for frame in traj:
        coords = frame.coordinates()
        r1 = np.linalg.norm(coords[1] - coords[7])
        r2 = np.linalg.norm(coords[7] - coords[2])

        dists.append(r1 - r2)

    plt.plot(dists, 627.509 * (np.array(energies) - min(energies)),
             lw=1.5,
             marker='o',
             label=label,
             color=color)

    return None


if __name__ == '__main__':

    fixed_k_energies = [-118.238095883376, -118.23552394357, -118.205142117475, -118.219337702456, -118.242397289654, -118.244974238858]
    adapt_k_energies = [-118.238095883376, -118.222903772576, -118.180220827535, -118.184871802568, -118.217284344293, -118.244974238858]
    irc_energies = [-118.244043, -118.244033, -118.243985, -118.243924, -118.243848, -118.243757, -118.243658, -118.243619, -118.243530, -118.242819, -118.241429, -118.239128, -118.235827, -118.231223, -118.224805, -118.216761, -118.207463, -118.197616, -118.187753, -118.177600, -118.176840, -118.177625, -118.187743, -118.197724, -118.207577, -118.216592, -118.223863, -118.228932, -118.232189, -118.234500, -118.236130, -118.237208, -118.237463, -118.237527, -118.237662, -118.237788, -118.237886, -118.237964, -118.237993, -118.238000, -118.238012]

    plot_line(fixed_k_energies,
              traj=gt.Data('neb_optimised_fixed.xyz'),
              label='fixed $k$',
              color=blues(0.5))

    plot_line(adapt_k_energies,
              traj=gt.Data('neb_optimised_adaptive.xyz'),
              label='adaptive $k$',
              color=blues(0.8))

    plot_line(irc_energies,
              traj=gt.Data('irc_IRC_Full.trj.xyz'),
              label='IRC',
              color='k')

    plt.legend()
    plt.ylabel('∆$E$ / kcal mol$^{-1}$')
    plt.xlabel('$r_1 - r_2$')
    plt.tight_layout()
    plt.savefig('tmp.pdf')
