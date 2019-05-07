#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "tatianefx"

from modules.chemistry.bond import Bond
from modules.chemistry.atom import Atom
from modules.chemistry.amino_acid import AminoAcid
from modules.common.helper import translation, resulting_vector_with_distance
from modules.common.config import PEPTIDE_BOND_LENGTH


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
                amino_acid = Molecule.__update_position(last_amino_acid.atoms[-1], amino_acid)
            peptide.append(amino_acid)

        return peptide

    @staticmethod
    def __update_position(atom_a: Atom, amino_acid: AminoAcid):
        matrix = []
        for atom in amino_acid.atoms:
            matrix.append(atom.position)

        v = resulting_vector_with_distance(atom_a.position, matrix[0], PEPTIDE_BOND_LENGTH)

        result = translation(matrix, v)

        # update positions
        for i in range(len(amino_acid.atoms)):
            amino_acid.atoms[i].position = result[i]

        return amino_acid


