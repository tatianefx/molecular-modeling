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
    print("\nThe best geometry:\n")
    print(genetic.the_best_individual.geometry)


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

# Gly-Gly-Gly   Total Energy = -696.5010784885544126
#
# geometry = "O -0.5150  1.1856  1.3202\n" \
# "O  1.9894 -1.0724 -0.8496\n" \
# "O -2.1156 -1.5861 -0.2908\n" \
# "O -4.2185 -1.2858  0.5064\n" \
# "N  1.8723  1.0949 -0.0233\n" \
# "N -1.6468  0.9523 -0.6943\n" \
# "N  4.5348 -1.3138  0.1518\n" \
# "C  0.6922  1.5116 -0.7467\n" \
# "C -0.5295  1.1952  0.0913\n" \
# "C -2.9510  0.6534 -0.1472\n" \
# "C  2.4263 -0.1719 -0.1363\n" \
# "C  3.6374 -0.3407  0.7580\n" \
# "C -3.1761 -0.8222  0.0606\n" \
# "H  0.7472  2.5920 -0.9060\n" \
# "H  0.6550  0.9970 -1.7116\n" \
# "H  2.2885  1.7640  0.6178\n" \
# "H -1.5581  0.9766 -1.7060\n" \
# "H -3.6937  1.0238 -0.8593\n" \
# "H -3.0622  1.1616  0.8151\n" \
# "H  3.3036 -0.6891  1.7400\n" \
# "H  4.1647  0.6096  0.8860\n" \
# "H  4.8156 -0.9965 -0.7752\n" \
# "H  4.0483 -2.1996  0.0196\n" \
# "H -2.2864 -2.5410 -0.1453\n" \
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
