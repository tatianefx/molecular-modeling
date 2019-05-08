#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "tatianefx"

from modules.chemistry.molecule import Molecule
from modules.chemistry.amino_acid import AminoAcid, Atom
from modules.genetic_algorithms.genetic_algorithm import GeneticAlgorithm
import sys
import psi4


amino_acid = input("Digite uma sequência de aminoácidos: ")
if AminoAcid.is_valid(amino_acid):
    peptide = Molecule.generate_peptide_molecule(amino_acid)
    genetic = GeneticAlgorithm(peptide)
    print("\n\nThe best fitness: " + str(genetic.the_best_individual.fitness))


# Ala-Gly     Total Energy = -528.7141214969851717
#
# geometry = "O  0.7678  1.0333  1.1238\n" \
# "O -1.8033 -1.3184  0.3291\n" \
# "O -3.5624  0.1136  0.2770\n" \
# "N -0.1435  0.2970 -0.8788\n" \
# "N  3.2928  0.5781  0.1123\n" \
# "C  2.1762 -0.1982 -0.4209\n" \
# "C  0.8780  0.4453  0.0498\n" \
# "C  2.2573 -1.6402  0.0644\n" \
# "C -1.4706  0.8298 -0.6749\n" \
# "C -2.3923 -0.1404  0.0183\n" \
# "H  2.2275 -0.1701 -1.5150\n" \
# "H  1.4277 -2.2305 -0.3403\n" \
# "H  2.1960 -1.7024  1.1571\n" \
# "H  3.1926 -2.1132 -0.2538\n" \
# "H  0.0278 -0.2295 -1.7304\n" \
# "H -1.8824  1.0584 -1.6618\n" \
# "H -1.4058  1.7404 -0.0721\n" \
# "H  4.1724  0.1652 -0.1968\n" \
# "H  3.3017  0.5174  1.1303\n" \
# "H -2.4226 -1.9313  0.7797\n"
#
# try:
#     psi4.core.set_output_file('output.dat', False)
#     psi4.set_memory('500 MB')
#
#     psi4.geometry(geometry)
#
#     psi4.energy('scf/cc-pvdz')
#     total_energy = psi4.core.get_variable('SCF TOTAL ENERGY')
# except RuntimeError:
#     pass

print("THE END")
sys.exit(0)
