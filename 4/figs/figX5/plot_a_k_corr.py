import matplotlib.pyplot as plt
import numpy as np
plt.style.use('paper')


if __name__ == '__main__':

    ks = [1.278, 0.718, 0.911, 1.386, 0.655, 0.707, 1.692, 0.796, 0.524]
    as_ = [2.705, 2.500, 2.449, 2.832, 2.562, 2.809, 2.688, 2.473, 2.542]

    print(np.average(ks), np.std(ks))
    print(np.average(as_), np.std(as_))


    plt.scatter(ks, as_)
    plt.scatter(ks[::3], as_[::3])
    plt.scatter(ks[1::3], as_[1::3])
    plt.scatter(ks[2::3], as_[2::3])

    plt.savefig('tmp.pdf')
