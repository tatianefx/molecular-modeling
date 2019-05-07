#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "tatianefx"

from modules.chemistry.molecule import Molecule
from modules.chemistry.amino_acid import AminoAcid, Atom
from modules.genetic_algorithms.genetic_algorithm import GeneticAlgorithm
import sys

amino_acid = input("Digite uma sequência de aminoácidos: ")
if AminoAcid.is_valid(amino_acid):
    peptide = Molecule.generate_peptide_molecule(amino_acid)
    genetic = GeneticAlgorithm(peptide)
    print(genetic.the_best_individual.fitness)

print("THE END")
sys.exit(0)
