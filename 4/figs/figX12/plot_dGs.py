import matplotlib.pyplot as plt
plt.style.use('paper')


if __name__ == '__main__':

    pbe0_d3bj_gp_svp = [-12.074774311498565, -10.89787291968087]

    # calc_dGcont(method='grimme')
    dG_conts_IGM = [13.544047928615218, 11.529058315582295]

    # calc_dGcont(method='ew')  # includes Grimme vib entropy
    dG_conts_EW_qrrho = [8.070701248421408, 6.055711635388487]

    # calc_dGcont(method='ew_rrho')  # includes HO vib entropy
    dG_conts_EW_rrho = [6.237131770347517, 2.3092743352067813]

    for i, dg_conts in enumerate([dG_conts_IGM,
                                  dG_conts_EW_rrho,
                                  dG_conts_EW_qrrho]):

        dG_N = pbe0_d3bj_gp_svp[0] + dg_conts[0]
        dG_CH =pbe0_d3bj_gp_svp[1] + dg_conts[1]

        plt.bar(x=[i*1.2 + 0.5],
                width=0.2,
                height=[dG_N],
                linewidth=0.8,
                label='IGM',
                color='tab:blue',
                alpha=1-(3-i)/10)

        plt.bar(x=[i*1.2 + 0.7],
                width=0.2,
                height=[dG_CH],
                linewidth=0.8,
                label='IGM',
                color='tab:orange',
                alpha=1-(3-i)/10)

    plt.plot([0, 10], [0, 0], ls='-', c='k', lw=0.8)

    plt.plot([0, 10], [-4.1, -4.1], ls='--', c='tab:blue')
    plt.plot([0, 10], [-5.3, -5.3], ls='--', c='tab:orange')

    plt.xticks([])
    plt.ylabel('$\Delta G_{calc}$ / kcal mol$^{-1}$')
    plt.xlim(0, 3.6)

    plt.tight_layout()
    plt.savefig('dgs.pdf')



