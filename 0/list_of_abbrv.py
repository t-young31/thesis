"""
Generate an alphabetically ordered list of abbreviations used, as a .tex file
"""
glossary = {
        'WF': 'Wave function',
        'HF': 'Hartree-Fock',
        'DFT': 'Density Functional Theory',
        'KS' : 'Kohn-Sham',
        'PT': 'Perturbation theory',
        'MP2': 'Moller-Plesset perturbation theory',
        'CC': 'Coupled cluster',
        'CI': 'Configuration interaction',
        'CC': 'Coupled cluster',
        'CCSD': 'Coupled cluster singles and doubles',
        'CCSD(T)': 'CCSD with perturbative triples',
        'MO': 'Molecular orbtial',
        'AO': 'Atomic orbtial',
        'GTO': 'Gaussian-type orbtial',
        'STO': 'Slater-type orbtial',
        'BSIE': 'Basis set incompleteness error',
        'BSSE': 'Basis set superposition error',
        'RI': 'Resolution of the identity',
        'LCC': 'Local coupled cluster',
        'DLPNO': 'Domain-based local pair natural orbital',
        'PP': 'Pseudopotential',
        'ECP': 'Effective core potential',
        'PES': 'Potential energy surface',
        'FEP': 'Free energy perturbation',
        'MD': 'Molecular dynamics',
        'MC': 'Metropolis Monte Carlo',
        'AIMD': 'Ab initio molecular dynamics',
        'ML': 'Machine learning',
        'QSAR': 'Quantitative structure activity relationships',
        'QM': 'Quantum mechanics',
        'RMSD': 'Root mean squared deviation',
        'NCI': 'Non-covalent interaction',
        'MOF': 'Metal organic framework',
        'API': 'Application programming interface',
        'GUI': 'Graphical user interface',
        'SSD': 'Sum squared distance',
        'ETKDG': 'Experimental-Torsion Distance Geometry with Knowlege',
        'ETDG': 'Experimental-Torsion Distance Geometry',
        'BFGS': 'Broyden-Fletcher-Goldfarb-Shanno optimisation',
        'ESP': 'Electrostatic potential',
        'COM': 'Centre of mass',
        'FF': 'Force field',
        'SMILES': 'Simplified molecular-input line-entry system',
        'MSD': 'Mean signed deviation',
        'HVTS': 'High-throughput virtual screening',
        'TS': 'Transition state',
        'AINR': 'Ab initio nanoreactor',
        'AFIR': 'Artificial force induced reaction',
        'GSM': 'Growing string method',
        'FSM': 'Freezing string method',
        'TSA': 'Transition state analogue',
        'NEB': 'Nudged elastic band',
        'MEP': 'Minimum energy pathway',
        'LEPS': 'London-Eyring-Polanyi-Sato',
        'NN': '[Artificial] neural network',
        'ANN': 'Artificial neural network',
        'ACSF': 'Atomic-centred symmetry functions',
        'GPR': 'Gaussian process regression',
        'MLP': 'Machine learned potential',
        'QSAR': 'Quantitative structure activity relationships',
        'HOMO': 'Highest occupied molecular orbtial',
        'LUMO': 'Lowest unoccupied molecular orbtial',
        'GAP': 'Gaussian approximation potential',
        'AIMD': 'Ab initio molecular dynamics',
        'PES': 'Potential energy surface',
        'MC': '(Metropolis) Monte Carlo',
        'GDML': 'Gradient-domain machine learning',
        'AL': 'Active learning',
        'RMSE': 'Root mean squared error',
        'MSE': 'Mean squared error',
        'MAD': 'Mean absolute deviation',
        'RDF': 'Radial distribution function',
        'VRI': 'Valley-ridge inflection point',
        'QM/MM': 'Quantum mechanics in a molecular mechanics environment',
        'AE': 'Absolute error',
        'NBO': 'Natural bond order',
        'DA': 'Diels-Alder',
        'EDDM': 'Electron density difference map',
        'MOF': 'Metal-organic frameworks',
        'PF': 'Partition function',
        'IGM': 'Ideal gas method',
        'PIB': 'Particle in a box',
        'RR': 'Ridgid rotor',
        'HO': 'Harmonic oscillator',
        'LFM': 'Low frequency mode',
        'DCM': 'Dichloromethane',
        'RESP': 'Restrained electrostatic potential',
        'ZPE': 'Zero point energy',
        'LCC': 'Local coupled cluster'
    }


def print_tex():
    """Print the file"""

    with open('list_of_abbrv.tex', 'w') as out:

        print(r'\documentclass[../main.tex]{subfiles}',
              r'\begin{document}',
              r'\begin{center}',
              r'		{\bfseries\Large \textsf{List of Abbreviations}}',
              r'\end{center}',
              r'\addcontentsline{toc}{chapter}{List of Abbreviations}',
              r'\begin{table}[h!]',
              r'\def\arraystretch{2.0}',
              r'\begin{tabularx}{\textwidth}{YY}',
              #r'\toprule',
              r'Abbreviation & Meaning \\',
              r'\hline',
              sep='\n', file=out)

        for i, (abbr, meaning) in enumerate(sorted(glossary.items(),
                                                   key=lambda x: x[0])):

            if i > 0 and i % 19 == 0:
                print(r'\end{tabularx}',
                      r'\end{table}',
                      r'\newpage',
                      r'\begin{table}[h!]',
                      r'\def\arraystretch{2.0}',
                      r'\begin{tabularx}{\textwidth}{YY}',
                      r'Abbreviation & Meaning \\',
                      r'\hline',
                      sep='\n', file=out)

            print(abbr + '\t&\t' + meaning + '\t' + r'\\', file=out)

        print(#r'\bottomrule',
              r'\end{tabularx}',
              r'\end{table}',
              r'\clearpage',
              r'\end{document}',
              sep='\n', file=out)

    return


if __name__ == '__main__':

    print_tex()
