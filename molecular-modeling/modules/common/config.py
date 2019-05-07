#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "tatianefx"

import os

dirname = os.path.dirname(__file__)
DIRECTORY_PATH = os.path.join(dirname, 'amino_acids/')

"""
'Typically the peptide bond(C-N) has a length of 1.32 to 1.33 Å.'

Advanced Biotechnology by R C Dubey
https://books.google.com.br/books?id=SKgrDAAAQBAJ&pg=PA457&lpg=PA457&dq=peptide+bond+is+typically+1.32&source=bl&ots=FPMGytIfNn&sig=ACfU3U3_949XhtXwkLUurCSbD4kKWGozKw&hl=pt-BR&sa=X&ved=2ahUKEwjD6bH88oniAhUFGbkGHbdHAoMQ6AEwC3oECAkQAQ#v=onepage&q=peptide%20bond%20is%20typically%201.32&f=false
"""
# BOND DISTANCE IN ANGSTRON
PEPTIDE_BOND_LENGTH = 1.32

# GENETIC ALGORITHM PARAMETERS
POPULATION_SIZE = 100
GENERATIONS = 100
MUTATION_RATE = 0.2
CROSSOVER_RATE = 0.2

