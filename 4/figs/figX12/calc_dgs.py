from otherm import Molecule, Constants


def calc_dEs():
    # SP M062X def2-TZVP RIJCOSX def2/J CPCM SMDsolvent "CH2Cl2"
    de_n_tz = -4215.709782735769 - (-3834.245108141723 + -381.459777502851)
    de_ch_tz = -4151.559879877889 - (-3770.093879006371 + -381.459777502851)

    # SP M062X def2-TZVP RIJCOSX def2/J CPCM SMDsolvent "CH2Cl2" GCP(DFT/TZ)
    de_n_tz_gp = -4215.634849988669 - (-3834.177956103723 + -381.454331078551)
    de_ch_tz_gp = -4151.483450930033 - (-3770.02499041277 + -381.454331078551)

    # SP M062X def2-SVP RIJCOSX def2/J CPCM SMDsolvent "CH2Cl2" GCP(DFT/SVZ)
    de_n_svp_gp = -4210.743598664712 - (-3829.760606254353 + -380.977592620334)
    de_ch_svp_gp = -4146.675153249457 - (
                -3765.689165730293 + -380.977592620334)

    # SP PBE0 D3BJ def2-SVP RIJCOSX GCP(DFT/SVP) def2/J CPCM
    de_n_svp_pbe0_gp = -4208.236872527182 - (-3827.494536207277 + -380.723093927666)
    de_ch_svp_pbe0_gp = -4144.216228540279 - (-3763.475767733541 + -380.723093927666)

    for (de_n, de_ch) in [(de_n_tz, de_ch_tz),
                          (de_n_tz_gp, de_ch_tz_gp),
                          (de_n_svp_gp, de_ch_svp_gp),
                          (de_n_svp_pbe0_gp, de_ch_svp_pbe0_gp)]:
        print('∆E_N', 627.509 * de_n)
        print('∆E_CH', 627.509 * de_ch)

    return


def calc_dGcont(method):

    ch_quinone = Molecule('CHcage_quinone_tmp.out', is_ts=False, real_freqs=True)
    ch = Molecule('CHcage_sp.out', is_ts=False, real_freqs=True)
    n_quinone = Molecule('Ncage_quinone_sp.out', is_ts=False, real_freqs=True)
    n = Molecule('Ncage_sp.out', is_ts=False, real_freqs=True)
    quinone = Molecule('quinone.out', is_ts=False,  real_freqs=True)

    for mol in (ch_quinone, ch, n_quinone, n, quinone):
        mol.calculate_thermochemistry(temp=298.15,
                                      ss='1M',
                                      method=method,
                                      w0=100,
                                      alpha=4,
                                      calc_sym=mol.n_atoms < 20)

    print('∆G_cont^IGM(Ncage)',
          Constants.j_to_kcal * (n_quinone.g_cont - (n.g_cont + quinone.g_cont)))

    print('∆G_cont^IGM(CHcage)',
          Constants.j_to_kcal * (ch_quinone.g_cont - (ch.g_cont + quinone.g_cont)))
    return


if __name__ == '__main__':

    # calc_dEs()
    # dEs = [-3.0729687985794034, -3.9052198486873904, -1.608184078096546,
    #        -2.591259956709946, -3.3884168387769753, -5.267874570149527]

    m062x_solv_gp_svp = [-3.3884168387769753, -5.267874570149527]
    pbe0_d3bj_gp_svp = [-12.074774311498565, -10.89787291968087]

    # calc_dGcont(method='grimme')
    dG_conts_IGM = [13.544047928615218, 11.529058315582295]

    # calc_dGcont(method='ew')  # includes Grimme vib entropy
    dG_conts_EW_qrrho = [8.070701248421408, 6.055711635388487]

    # calc_dGcont(method='ew_rrho')  # includes HO vib entropy
    dG_conts_EW_rrho = [6.237131770347517, 2.3092743352067813]

    for dg_conts in (dG_conts_IGM, dG_conts_EW_rrho, dG_conts_EW_qrrho):
        print('∆G_N = ', pbe0_d3bj_gp_svp[0] + dg_conts[0])
        print('∆G_CH = ', pbe0_d3bj_gp_svp[1] + dg_conts[1])
