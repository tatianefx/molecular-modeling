#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "tatianefx"

from enum import Enum
from modules.common.helper import calculates_distance_a_b
from modules.chemistry.atom import Atom


class BondType(Enum):
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    PEPTIDE = 4


class Bond:

    def __init__(self, first_atom: Atom, second_atom: Atom, bond_type):
        self.first_atom = first_atom
        self.second_atom = second_atom
        self.bond_type = bond_type
        self.bond_length = calculates_distance_a_b(first_atom.position, second_atom.position)

