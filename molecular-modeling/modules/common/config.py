#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "tatianefx"

import os

dirname = os.path.dirname(__file__)
DIRECTORY_PATH = os.path.join(dirname, 'amino_acids/')

# BOND DISTANCE IN ANGSTRON
PEPTIDE_BOND_LENGTH = 1

# GENETIC ALGORITHM PARAMETERS
POPULATION_SIZE = 100
GENERATIONS = 100
MUTATION_RATE = 0.2
CROSSOVER_RATE = 0.2

