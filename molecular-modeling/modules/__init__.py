from modules.chemistry.molecule import Molecule
from modules.chemistry.amino_acid import AminoAcid, Atom
import psi4

amino_acid = input("Digite uma sequência de aminoácidos: ")
geometry = ""

if AminoAcid.is_valid(amino_acid):
    peptide = Molecule.generate_peptide_molecule(amino_acid)
    for item in peptide:
        for atom in item.atoms:
            geometry = geometry + atom.acronym.value + ' '
            geometry = geometry + str(atom.position[0]) + ' '
            geometry = geometry + str(atom.position[1]) + ' '
            geometry = geometry + str(atom.position[2]) + '\n'
    print(geometry)


psi4.core.set_output_file('output.dat', False)
psi4.set_memory('500 MB')

h2o = psi4.geometry(geometry)

psi4.energy('scf/cc-pvdz')

