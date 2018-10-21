from enum import Enum
from modules.chemistry.atom import Atom, AtomName
from modules.chemistry.bond import Bond, BondType
from modules.common.config import DIRECTORY_PATH
import random


ALPHABET = "ARNDCEQGHILKMFPSTWYVUOX"


class AminoAcidType(Enum):
    ALANINE = ('A', DIRECTORY_PATH + "Alanine.mol")
    ARGININE = ('R', DIRECTORY_PATH + "Arginine.mol")
    ASPARAGINE = ('N', DIRECTORY_PATH + "Asparagine.mol")
    ASPARTIC_ACID = ('D', DIRECTORY_PATH + "Aspartic_acid.mol")
    CYSTEINE = ('C', DIRECTORY_PATH + "Cysteine.mol")
    GLUTAMIC_ACID = ('E', DIRECTORY_PATH + "Glutamic_acid.mol")
    GLUTAMINE = ('Q', DIRECTORY_PATH + "Glutamine.mol")
    GLYCINE = ('G', DIRECTORY_PATH + "Glycine.mol")
    HISTIDINE = ('H', DIRECTORY_PATH + "Histidine.mol")
    ISOLEUCINE = ('I', DIRECTORY_PATH + "Isoleucine.mol")
    LEUCINE = ('L', DIRECTORY_PATH + "Leucine.mol")
    LYSINE = ('K', DIRECTORY_PATH + "Lysine.mol")
    METHIONINE = ('M', DIRECTORY_PATH + "Methionine.mol")
    PHENYLALANINE = ('F', DIRECTORY_PATH + "L-phenylalanine.mol")
    PROLINE = ('P', DIRECTORY_PATH + "Proline.mol")
    SERINE = ('S', DIRECTORY_PATH + "Serine.mol")
    THREONINE = ('T', DIRECTORY_PATH + "Threonine.mol")
    TRYPTOPHAN = ('W', DIRECTORY_PATH + "Tryptophan.mol")
    TYROSINE = ('Y', DIRECTORY_PATH + "Tyrosine.mol")
    VALINE = ('V', DIRECTORY_PATH + "Valine.mol")
    SELENOCYSTEINE = ('U', DIRECTORY_PATH + "Selenocysteine.mol")
    PYRROLYSINE = ('O', DIRECTORY_PATH + "Pyrrolysine.mol")
    UNKNOWN = ('X', '')

    @staticmethod
    def list():
        return list(map(lambda a: a.value, AminoAcidType))

    @staticmethod
    def random_choice():
        values = AminoAcidType.list()
        values.pop()
        return AminoAcidType(random.choice(values))


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
        if amino_acid == AminoAcidType.ALANINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.ALANINE)
        elif amino_acid == AminoAcidType.ARGININE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.ARGININE)
        elif amino_acid == AminoAcidType.ASPARAGINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.ASPARAGINE)
        elif amino_acid == AminoAcidType.ASPARTIC_ACID.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.ASPARTIC_ACID)
        elif amino_acid == AminoAcidType.CYSTEINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.CYSTEINE)
        elif amino_acid == AminoAcidType.GLUTAMIC_ACID.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.GLUTAMIC_ACID)
        elif amino_acid == AminoAcidType.GLUTAMINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.GLUTAMINE)
        elif amino_acid == AminoAcidType.GLYCINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.GLYCINE)
        elif amino_acid == AminoAcidType.HISTIDINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.HISTIDINE)
        elif amino_acid == AminoAcidType.ISOLEUCINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.ISOLEUCINE)
        elif amino_acid == AminoAcidType.LEUCINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.LEUCINE)
        elif amino_acid == AminoAcidType.LYSINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.LYSINE)
        elif amino_acid == AminoAcidType.METHIONINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.METHIONINE)
        elif amino_acid == AminoAcidType.PHENYLALANINE.value:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.PHENYLALANINE)
        elif amino_acid == AminoAcidType.PROLINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.PROLINE)
        elif amino_acid == AminoAcidType.SERINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.SERINE)
        elif amino_acid == AminoAcidType.THREONINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.THREONINE)
        elif amino_acid == AminoAcidType.TRYPTOPHAN.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.TRYPTOPHAN)
        elif amino_acid == AminoAcidType.TYROSINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.TYROSINE)
        elif amino_acid == AminoAcidType.VALINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.VALINE)
        elif amino_acid == AminoAcidType.SELENOCYSTEINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.SELENOCYSTEINE)
        elif amino_acid == AminoAcidType.PYRROLYSINE.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.PYRROLYSINE)
        elif amino_acid == AminoAcidType.UNKNOWN.value[0]:
            return AminoAcid.__get_amino_acid_by_path(AminoAcidType.random_choice())
        else:
            return

    @staticmethod
    def __get_amino_acid_by_path(type: AminoAcidType):
        with open(type.value[1]) as f:
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

