from modules.chemistry.molecule import Molecule
from modules.chemistry.amino_acid import AminoAcid


amino_acid = input("Digite uma sequência de aminoácidos: ")

if AminoAcid.is_valid(amino_acid):
    Molecule.generate_peptide_molecule(amino_acid)
