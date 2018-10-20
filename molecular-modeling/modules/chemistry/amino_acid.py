from enum import Enum
from modules.chemistry.atom import Atom
from modules.chemistry.bond import Bond
from modules.chemistry.molecule import Molecule


ALPHABET =  "ARNDCEQGHILKMFPSTWYVUOX"
DIRECTORY_PATH = "/Users/tatianefx/PycharmProjects/molecular-modeling/modules/common/amino_acids/"


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
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Alanine.mol")
        elif amino_acid == AminoAcidType.ARGININE.value:
            return AminoAcid.__get_amino_acid_by_path(DIRECTORY_PATH + "Alanine.mol")
        return

    @staticmethod
    def __get_amino_acid_by_path(amino_acid_path):
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
                position = []
                position.append(float(array[0]))
                position.append(float(array[1]))
                position.append(float(array[2]))
                acronym = array[3]
                atoms.append(Atom(acronym, position))

            #  Read bonds
            for i in range(num_bonds):
                array = next(f).split()
                first_atom = int(array[0])
                second_atom = int(array[1])
                bond_type = int(array[2])
                bonds.append(Bond(atoms[first_atom-1], atoms[second_atom-1], bond_type))

            return Molecule(atoms, bonds)
        return None

