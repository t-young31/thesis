from cgbind import Linker, Cage, Config, XTB
from cgbind.input_output import xyzs2xyzfile, xyzfile2xyzs
import data
import os
from subprocess import Popen, PIPE

Config.n_cores = 8

if __name__ == '__main__':

    for struct in data.structs:
        basename = struct.stuct_filename.replace('.mol2', '')

        crystal_strcut_mol2_filename = f'{basename}.mol2'
        crystal_strcut_xyz_filename = f'{basename}.xyz'

        # Convert .mol2 to .xyz with open babel so RMSD can be used
        sp = Popen(['obabel', '-imol2', crystal_strcut_mol2_filename, '-oxyz', '-O', crystal_strcut_xyz_filename],
                   stdout=PIPE, stderr=PIPE)
        sp.wait()

        if True:
            # Construct the cage using cgbind
            if basename == 'RIJGUM_edit':
                linker = Linker(arch_name=struct.arch, xyzs=xyzfile2xyzs('fujita_linker.xyz'))

            else:
                linker = Linker(arch_name=struct.arch, smiles=struct.linker_smlies, n_confs=1000, use_etdg_confs=False)

            xyzs2xyzfile(xyzs=linker.xyzs, basename=f'{basename}_linker')
            cage = Cage(linker=linker, metal=struct.metal, metal_charge=struct.metal_charge)
            cage.optimise(method=XTB, n_cores=8)
            os.remove('cage_linker_opt_xtb.out')
            xyzs2xyzfile(xyzs=cage.xyzs, basename=f'{basename}_cgbind')

        aligned_crystal_strcut_xyz_filename = crystal_strcut_xyz_filename.replace('.xyz', '_aligned.xyz')
        if os.path.exists(aligned_crystal_strcut_xyz_filename):
            sp = Popen(['calculate_rmsd', '--reorder', aligned_crystal_strcut_xyz_filename, f'{basename}_cgbind.xyz'], stdout=PIPE, stderr=PIPE)
        else:
            sp = Popen(['calculate_rmsd', '--reorder', crystal_strcut_xyz_filename, f'{basename}_cgbind.xyz'], stdout=PIPE, stderr=PIPE)
        stdout, stderr = sp.communicate()
        rmsd = stdout.decode("utf-8").split()[0]
        print(rmsd)
