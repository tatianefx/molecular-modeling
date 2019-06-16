#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "tatianefx"

import os

dirname = os.path.dirname(__file__)
DIRECTORY_PATH = os.path.join(dirname, 'amino_acids/')

"""
'It is clear that the bond length for the peptide bond(1.32Å) is shorter than that for the single N-C bond(1.45Å).

Exploring Protein Structure: Principles and Practice
by Tim Skern
https://books.google.com.br/books?id=lQtjDwAAQBAJ&pg=PA63&dq=peptide+bond+length+1.32&hl=pt-BR&sa=X&ved=0ahUKEwiL1-zDioziAhUlLLkGHdIjDTMQ6AEIKTAA#v=onepage&q=peptide%20bond%20length%201.32&f=false
"""
# BOND DISTANCE IN ANGSTRON
PEPTIDE_BOND_LENGTH = 1.32

# GENETIC ALGORITHM PARAMETERS
POPULATION_SIZE = 5
GENERATIONS = 5
MUTATION_RATE = 0.2
CROSSOVER_RATE = 0.2

