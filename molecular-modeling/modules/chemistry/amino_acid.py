from enum import Enum
from modules.chemistry.atom import Atom, AtomName
from modules.chemistry.bond import Bond, BondType


ALPHABET =  "ARNDCEQGHILKMFPSTWYVUOX"
DIRECTORY_PATH = "/Users/tatianefx/PycharmProjects/tcc-molecular-modeling/molecular-modeling/modules/common/amino_acids/"


class AminoAcidType(Enum):
    ALANINE = 'A'
    ARGININE = 'R'
    ASPARAGINE = 'N'
    ASPARTIC_ACID = 'D'
    CYSTEINE = 'C'
    GLUTAMIC_ACID = 'E'
    GLUTAMINE = 'Q'
    GLYCINE = 'G'
    HISTIDINE = 'H'
    ISOLEUCINE = 'I'
    LEUCINE = 'L'
    LYSINE = 'K'
    METHIONINE = 'M'
    PHENYLALANINE = 'F'
    PROLINE = 'P'
    SERINE = 'S'
    THREONINE = 'T'
    TRYPTOPHAN = 'W'
    TYROSINE = 'Y'
    VALINE = 'V'
    SELENOCYSTEINE = 'U'
    PYRROLYSINE = 'O'
    UNKNOWN = 'X'


class AminoAcid:

    def __init__(self, atoms: [Atom], bonds: [Bond], amino_acid_type):
        self.atoms = atoms
        self.bonds = bonds
        self.amino_acid_type = amino_acid_type

    def eliminate_hydroxyl_radical(self):
        bond = self.bonds[0]
        carbon = bond.second_atom
        del self.bonds[0]
        del self.bonds[0]
        return carbon

    def eliminate_hydrogen(self, carbon):
        for (i, bond) in enumerate(self.bonds):
            first_atom = bond.first_atom
            second_atom = bond.second_atom
            if first_atom.acronym == AtomName.N and second_atom.acronym == AtomName.H:
                bond.second_atom = carbon
                bond.bond_type = BondType.PEPTIDE
                self.bonds[i] = bond
                break


    @staticmethod
    def is_valid(amino_acid_sequence):
        for c in amino_acid_sequence:
            if c not in ALPHABET:
                print("Aminoácido desconhecido: " + c)
                return False
        return True

    @staticmethod
    def generate_amino_acid_molecule(amino_acid):
        if amino_acid == AminoAcidType.ALANINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Alanine.mol", AminoAcidType.ALANINE)
        elif amino_acid == AminoAcidType.ARGININE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Arginine.mol", AminoAcidType.ARGININE)
        elif amino_acid == AminoAcidType.ASPARAGINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Asparagine.mol", AminoAcidType.ASPARAGINE)
        elif amino_acid == AminoAcidType.ASPARTIC_ACID.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Aspartic_acid.mol", AminoAcidType.ASPARTIC_ACID)
        elif amino_acid == AminoAcidType.CYSTEINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Cysteine.mol", AminoAcidType.CYSTEINE)
        elif amino_acid == AminoAcidType.GLUTAMIC_ACID.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Glutamic_acid.mol", AminoAcidType.GLUTAMIC_ACID)
        elif amino_acid == AminoAcidType.GLUTAMINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Glutamine.mol", AminoAcidType.GLUTAMINE)
        elif amino_acid == AminoAcidType.GLYCINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Glycine.mol", AminoAcidType.GLYCINE)
        elif amino_acid == AminoAcidType.HISTIDINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Histidine.mol", AminoAcidType.HISTIDINE)
        elif amino_acid == AminoAcidType.ISOLEUCINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Isoleucine.mol", AminoAcidType.ISOLEUCINE)
        elif amino_acid == AminoAcidType.LEUCINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Leucine.mol", AminoAcidType.LEUCINE)
        elif amino_acid == AminoAcidType.LYSINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Lysine.mol", AminoAcidType.LYSINE)
        elif amino_acid == AminoAcidType.METHIONINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Methionine.mol", AminoAcidType.METHIONINE)
        elif amino_acid == AminoAcidType.PHENYLALANINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "L-phenylalanine.mol", AminoAcidType.PHENYLALANINE)
        elif amino_acid == AminoAcidType.SERINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Serine.mol", AminoAcidType.SERINE)
        elif amino_acid == AminoAcidType.THREONINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Threonine.mol", AminoAcidType.THREONINE)
        elif amino_acid == AminoAcidType.TRYPTOPHAN.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Tryptophan.mol", AminoAcidType.TRYPTOPHAN)
        elif amino_acid == AminoAcidType.TYROSINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Tyrosine.mol", AminoAcidType.TYROSINE)
        elif amino_acid == AminoAcidType.VALINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Valine.mol", AminoAcidType.VALINE)
        elif amino_acid == AminoAcidType.SELENOCYSTEINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Selenocysteine.mol", AminoAcidType.SELENOCYSTEINE)
        elif amino_acid == AminoAcidType.PYRROLYSINE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Pyrrolysine.mol", AminoAcidType.PYRROLYSINE)
        else:
            return

    @staticmethod
    def __get_amino_acid_by_path(amino_acid_path, type):
        with open(amino_acid_path) as f:
            size = 0
            num_bonds = 0

            atoms = []
            bonds = []

            # Read basic informations
            for (i, line) in enumerate(f):
                if i == 3:
                    array = line.split()
                    size = int(array[0])
                    num_bonds = int(array[1])
                    break

            #  Read atoms
            for i in range(size):
                array = next(f).split()
                position = list()
                position.append(float(array[0]))
                position.append(float(array[1]))
                position.append(float(array[2]))
                acronym = AtomName(array[3])
                atoms.append(Atom(acronym, position))

            #  Read bonds
            for i in range(num_bonds):
                array = next(f).split()
                first_atom = int(array[0])
                second_atom = int(array[1])
                bond_type = BondType(int(array[2]))
                bonds.append(Bond(atoms[first_atom-1], atoms[second_atom-1], bond_type))

            return AminoAcid(atoms, bonds, type)
        return None

