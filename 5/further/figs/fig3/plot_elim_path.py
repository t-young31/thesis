import matplotlib.pyplot as plt
from autode.constants import Constants
import numpy as np
plt.style.use('paper')


irc_energies = [-693.246013, -693.245961, -693.245904, -693.245880, -693.245778, -693.245477, -693.245108, -693.244620, -693.243974, -693.243260, -693.239943, -693.235541, -693.230949, -693.226579, -693.220614, -693.218199, -693.216409, -693.221722, -693.225620, -693.228950, -693.232605, -693.236353, -693.240124, -693.242712, -693.243990, -693.244901, -693.245051, -693.245112, -693.245273, -693.245443, -693.245603, -693.245746, -693.245871]
adapt_energies = [-693.203150916778, -693.18186905008, -693.175742753741, -693.175481244166, -693.192402185035, -693.232079484957, -693.242620582565]

if __name__ == '__main__':

    plt.plot(np.linspace(0, 6, num=len(irc_energies)),
             Constants.ha2kcalmol * (np.array(irc_energies) - min(irc_energies)),
             c='green')

    plt.plot(np.linspace(0, 6, num=len(adapt_energies)),
             Constants.ha2kcalmol * (np.array(adapt_energies) - min(irc_energies)),
             c='b')
    plt.savefig('elim_path.png', dpi=300)
