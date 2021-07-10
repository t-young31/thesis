import matplotlib.pyplot as plt
reds = plt.get_cmap('Reds')
blues = plt.get_cmap('Blues')
plt.style.use('paper')


if __name__ == '__main__':

    dg_igm, dg_stderr = (-8.641285419935713, 0.3977607904106817)

    ew_sigma = 27.319517253413203*298.15/ (1000 * 4.184)
    dg_ew_no_del, dg_ew_stderr_no_del = (0.3346700524075954, max(0.39781974977271223, ew_sigma))

    dg_ew, dg_ew_stderr = (0.3190447507612837, 0.3978197497727124)


    plt.figure(figsize=(5, 5))
    plt.bar(x=[0],
            width=0.5,
            height=[dg_igm],
            yerr=[dg_stderr],
            linewidth=0.8,
            edgecolor=reds(0.6),
            label='IGM',
            color=[reds(0.4)],
            capsize=2)

    # ∆G where rather than treating the translational part of one of the waters
    # in the octamer cluster with an EW, use the vibrational contribution and
    # do not delete freqencies
    plt.bar(x=[0.5],
            width=0.5,
            height=[dg_ew_no_del],
            yerr=[dg_ew_stderr_no_del],
            linewidth=0.8,
            edgecolor=blues(0.7),
            color=[blues(0.5)],
            label='exp$_1$',
            capsize=2)

    plt.bar(x=[1],
            width=0.5,
            height=[dg_ew],
            yerr=[dg_ew_stderr],
            linewidth=0.8,
            edgecolor=blues(0.9),
            color=[blues(0.7)],
            label='exp$_2$',
            capsize=2)

    plt.ylabel(r'$\Delta G$ / kcal mol$^{-1}$')
    plt.plot([-0.5, 4], [0, 0], c='k', lw=1.5)
    plt.xticks([])
    plt.legend()
    plt.xlim(-0.5, 1.5)
    plt.ylim(-12, 10)
    plt.tight_layout()
    plt.savefig('water_in_water_dGs.pdf')
