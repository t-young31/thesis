import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

# Just generated from cgbind
structs_and_rmsd = {
'EZEVAI': 0.4555875367756014,  # calculate_rmsd --reorder --no-hydrogen EZEVAI_edit_cgbind.xyz EZEVAI_edit.xyz
'OVUSIJ': 0.7954149721879067,  # calculate_rmsd --reorder --no-hydrogen OVUSIJ_edit_cgbind.xyz OVUSIJ_edit_aligned.xyz
'JIZPOZ': 0.6275429875993981,  # calculate_rmsd --reorder --no-hydrogen JIZPOZ_edit_cgbind.xyz JIZPOZ_edit_aligned.xyz
'GARWUR': 1.3783727565969857,  # calculate_rmsd --reorder --no-hydrogen GARWUR_edit_cgbind.xyz GARWUR_edit_aligned.xyz
'REDZEH': 0.8203224061058818,  # calculate_rmsd --reorder --no-hydrogen REDZEH_edit_cgbind.xyz REDZEH_edit_aligned.xyz
'SAYGUX': 1.099103627847991,   # calculate_rmsd --reorder  --no-hydrogen SAYGUX_edit_cgbind.xyz SAYGUX_edit_aligned.xyz
'GEGZAU': 1.2952912445009523,  # calculate_rmsd --reorder --no-hydrogen GEGZAU_edit_cgbind.xyz GEGZAU_edit.xyz
'RIJGUM': 1.2120886129389472,  # calculate_rmsd --reorder --no-hydrogen RIJGUM_edit_cgbind.xyz RIJGUM_edit_aligned.xyz
'OLOBUN': 1.4028090163046272   # calculate_rmsd --reorder --no-hydrogen OLOBUN_edit_cgbind.xyz OLOBUN_edit_aligned.xyz
}

# cgbind and XTB optimised
structs_and_rmsd_xtb = {
'EZEVAI': 0.4147376408064745,
'OVUSIJ': 0.9189640558267427,
'JIZPOZ': 0.4708358528345689,
'GARWUR': 2.0189383686743128,
'REDZEH': 0.2500992130492709,
'SAYGUX': 0.7356584908950148,
'GEGZAU': 1.735986490219473,
'RIJGUM': 0.4752107585490222,
'OLOBUN': 2.0563508944518905
}


# cgbind and ORCA optimised at: LooseOpt PBE RI D3BJ def2-SVP def2/J
structs_and_rmsd_orca = {
'EZEVAI': 0.46320740933695437,
'OVUSIJ': 1.7665076439055871,
'JIZPOZ': 0.4721954379311944,
'GARWUR': 2.8465234364911867,       # Not converged
'REDZEH': 0.6742699119961512,       # Not converged
'SAYGUX': 1.178039608916467,
'GEGZAU': 1.9847546985376538,       # Not converged
'RIJGUM': 0.6849640417339941,
'OLOBUN': 1.7655194847745486        # Not converged
}


if __name__ == '__main__':

    rmsd_matrix = np.zeros((3, 9))
    for i, val in enumerate(structs_and_rmsd.values()):
        rmsd_matrix[0, i] = val

    for i, val in enumerate(structs_and_rmsd_xtb.values()):
        rmsd_matrix[1, i] = val

    for i, val in enumerate(structs_and_rmsd_orca.values()):
        rmsd_matrix[2, i] = val

    print(f'MAD(cgbind) = {np.average(np.array(list(structs_and_rmsd.values()))):.2f} Å')
    print(f'MAD(XTB) = {np.average(np.array(list(structs_and_rmsd_xtb.values()))):.2f} Å')
    print(f'MAD(DFT) = {np.average(np.array(list(structs_and_rmsd_orca.values()))):.2f} Å')

    fig, ax = plt.subplots()
    heatplot = ax.imshow(rmsd_matrix, cmap='RdYlGn_r', vmin=0.0, vmax=2.0)
    ax.set_xticklabels([])

    ax.set_yticklabels([])
    ax.set_yticks([])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(heatplot, cax=cax)
    plt.tight_layout()
    # plt.show()
    plt.savefig('rmsd_fig.png', dpi=500)
