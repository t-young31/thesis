structs = []


class Struct:

    def __init__(self, arch, metal, metal_charge, linker_smlies, struct_filename):
        """
        Class to hold a binding affinity data-point.

        :param arch: (str) Metallocage architecture e.g. M2L4
        :param metal: (str) Label of the metal atom
        :param metal_charge: (int) Charge on the metal
        :param linker_smlies: (str) SMILES string of the cage linker
        :param substrate_smiles: (str) SMILES string of the substrate
        :param struct_filename: (str) Name of the .mol2 file containing the strucutre
        """

        self.arch = arch
        self.metal = metal
        self.metal_charge = metal_charge
        self.linker_smlies = linker_smlies
        self.stuct_filename = struct_filename


structs = [# Is template
            Struct(arch='m2l4', metal='Pd', metal_charge=2,
                   linker_smlies='C1(C#CC2=CC=CN=C2)=CC(C#CC3=CN=CC=C3)=CC=C1',
                   struct_filename='EZEVAI_edit.mol2'),
            Struct(arch='m2l4', metal='Pd', metal_charge=2,
                   linker_smlies='NC1=CC(C#CC2=CC(OC)=CN=C2)=CC(C#CC3=CN=CC(OC)=C3)=C1',
                   struct_filename='OVUSIJ_edit.mol2'),
            Struct(arch='m2l4', metal='Pd', metal_charge=2,
                   linker_smlies='COC1=CC(C#CC2=CC=CN=C2)=C(N)C(C#CC3=CN=CC=C3)=C1',
                   struct_filename='JIZPOZ_edit.mol2'),

            # Is template
            Struct(arch='m4l6', metal='Fe', metal_charge=4,
                   linker_smlies='O=C(NC1=CC=CC2=C1C=CC=C2NC(C3=C([O-])C([O-])=CC=C3)=O)C4=C([O-])C([O-])=CC=C4',
                   struct_filename='GARWUR_edit.mol2'),
            Struct(arch='m4l6n', metal='Mn', metal_charge=2,
                   linker_smlies='CC([C@H]1C2)(C)[C@@H](C1)C3=C2C=CC(C(C=C4)=NC=C4C5=CC=C(C6=CC=C(C[C@@H]7C(C)(C)[C@H]8C7)C8=N6)N=C5)=N3',
                   struct_filename='REDZEH_edit.mol2'),
            Struct(arch='m4l6n', metal='Fe', metal_charge=2,
                   linker_smlies='C1(/C=N/C2=CC=C(C=C2)C3=CC4=C(C=C3)C=C(C5=CC=C(/N=C/C6=NC=CC=C6)C=C5)C=C4)=NC=CC=C1',
                   struct_filename='SAYGUX_edit.mol2'),

            # Is template
            Struct(arch='m6l8', metal='Pd', metal_charge=2,
                   linker_smlies='O=C(NC1=CC=CN=C1)C2=CC(C(NC3=CC=CN=C3)=O)=CC(C(NC4=CC=CN=C4)=O)=C2',
                   struct_filename='GEGZAU_edit.mol2'),
            Struct(arch='m6l8', metal='Pd', metal_charge=2,
                   linker_smlies='C1(CC2=CC=NC=C2)=CC(CC3=CC=NC=C3)=CC(CC4=CC=NC=C4)=C1',
                   struct_filename='RIJGUM_edit.mol2'),
            Struct(arch='m6l8', metal='Cu', metal_charge=2,
                   linker_smlies='O=P(NC1=CN=CC=C1)(NC2=CC=CN=C2)NC3=CC=CN=C3',
                   struct_filename='OLOBUN_edit.mol2'),
            ]
