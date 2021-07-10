import gaptrain as gt
import numpy as np


if __name__ == '__main__':

    for frame in gt.Data('IRC_B.xyz'):
        coords = frame.coordinates()

        d_cf = np.linalg.norm(coords[0] - coords[2])
        d_ccl = np.linalg.norm(coords[1] - coords[2])

        print(round(d_ccl, 4))
