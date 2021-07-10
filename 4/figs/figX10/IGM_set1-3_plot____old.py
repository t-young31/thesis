import matplotlib.pyplot as plt
import numpy as np
plt.style.use('graph')


class oneATM_SS:

    plata2015_calc = [-171.8,
                      -218.2]

    valsov_calc = [-141.7,
                   -141.7,
                   -133.1,
                   -138.1,
                   -142.6,
                   -99.0,
                   -150.1,
                   -150.2,
                   -108.5,
                   -98.5
                   ]

    guthrie_calc = [-0.3,
                    -21.2,
                    5.6,
                    4.5,
                    3.7
                    ]



class oneM_SS:

    plata2015_calc = [
        -144.7,
        -190.8
    ]

    valsov_calc = [
        -113.3,
        -113.3,
        -105.6,
        -111.3,
        -115.8,
        -72.2,
        -123.4,
        -125.1,
        -81.8,
        -71.9
    ]

    guthrie_calc = [
        -0.4,
        -21.3,
        5.5,
        4.4,
        3.6
    ]


plata2015_expt = [-104.2,
                  -129.7]

valsov_expt = [-44,
               -69,
               -24.5,
               -80.8,
               -38.7,
               -28.4,
               -67.3,
               -51,
               -32.6,
               -18.4
               ]

guthrie_expt = [42.8,
                24.6,
                -9.7,
                16.9,
                46.3
                ]


def avg_mae(calc_list, expt_list):

    maes = []

    for i in range(len(calc_list)):
        for j in range(len(calc_list[i])):
            maes.append(np.abs(calc_list[i][j] - expt_list[i][j]))

    return np.average(np.array(maes))


if __name__ == '__main__':

    all_expt = plata2015_expt + valsov_expt + guthrie_expt
    min_val, max_val = -250, 50

    ss = oneM_SS

    print(avg_mae([ss.guthrie_calc, ss.plata2015_calc, ss.valsov_calc], [guthrie_expt, plata2015_expt, valsov_expt]))

    plt.scatter(guthrie_expt, ss.guthrie_calc, label='SET1', marker='+')
    plt.scatter(plata2015_expt, ss.plata2015_calc, label='SET2', marker='+')
    plt.scatter(valsov_expt, ss.valsov_calc, label='SET3', marker='+')
    plt.plot([min_val, max_val], [min_val, max_val], c='k', ls='--')
    plt.xlim(min_val, max_val)
    plt.ylim(min_val, max_val)
    plt.legend(prop={'weight':'bold'})
    plt.xlabel('$\Delta S_{expt}$ / J K$^{-1}$ mol$^{-1}$')
    plt.ylabel('$\Delta S_{calc}$ / J K$^{-1}$ mol$^{-1}$')
    plt.show()
