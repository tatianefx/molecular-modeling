from modules.chemistry.bond import Bond
from modules.chemistry.amino_acid import AminoAcid, AminoAcidType


class Molecule:

    def __init__(self, atoms: [AminoAcid], bonds: [Bond]):
        self.atoms = atoms
        self.bonds = bonds

    @staticmethod
    def generate_peptide_molecule(amino_acid_sequence):
        peptide = []

        for a in amino_acid_sequence:
            amino_acid = AminoAcid.generate_amino_acid_molecule(a)
            if len(peptide) > 0:
                last_amino_acid: AminoAcid = peptide[-1]
                carbon = last_amino_acid.eliminate_hydroxyl_radical()
                amino_acid.eliminate_hydrogen(carbon)
            peptide.append(amino_acid)

        return peptide

