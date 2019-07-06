#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "tatianefx"

from modules.chemistry.molecule import Molecule
from modules.chemistry.amino_acid import AminoAcid, Atom
from modules.genetic_algorithms.genetic_algorithm import GeneticAlgorithm
from modules.common.helper import normalize_vector
import sys
import psi4
import numpy as np


""" Molecular Modelling """

amino_acid = input("Enter an amino acid sequence: ")

if AminoAcid.is_valid(amino_acid):
    peptide = Molecule.generate_peptide_molecule(amino_acid)
    genetic = GeneticAlgorithm(peptide)
    print("\n\nThe best fitness: " + str(genetic.the_best_individual.fitness))
    print("\nThe best geometry:\n")
    print(genetic.the_best_individual.geometry)
    matrix = genetic.the_best_individual.matrix_2d()

# ============================================

""" Energy references """

# Ala-Gly     Total Energy = -528.7141214969851717

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

# Ala-Gly-Ser Total Energy = -849.4371377028904
#
# geometry =  "O -2.4488 -2.2878 -0.0662\n" \
#    "O  0.0028  1.2407 -1.2777\n" \
#    "O -4.8271  0.7290 -0.4805\n" \
#    "O  2.4196 -1.2046  0.7962\n" \
#    "O -3.9235 -0.0531  1.4491\n" \
#    "N -1.4203  0.3319  0.3192\n" \
#    "N  2.0891  1.0752  0.5210\n" \
#    "N  5.1365 -0.8175  0.5225\n" \
#    "C -2.5414  0.0608 -0.5520\n" \
#    "C  0.7346  1.1186  1.0215\n" \
#    "C  4.2279  0.1407 -0.1019\n" \
#    "C -2.4758 -1.3448 -1.1307\n" \
#    "C -0.2350  0.9029 -0.1208\n" \
#    "C  2.8278 -0.0982  0.4492\n" \
#    "C -3.8028  0.2403  0.2667\n" \
#    "C  4.2364 -0.0338 -1.6155\n" \
#    "H -2.5336  0.8015 -1.3597\n" \
#    "H  0.6090  0.3556  1.7960\n" \
#    "H  0.5598  2.1082  1.4534\n" \
#    "H -1.5100  0.0793  1.2996\n" \
#    "H  4.5602  1.1516  0.1582\n" \
#    "H -1.5649 -1.4826 -1.7219\n" \
#    "H -3.3392 -1.5633 -1.7676\n" \
#    "H  2.4969  1.9438  0.1878\n" \
#    "H  3.5634  0.6883 -2.0908\n" \
#    "H  3.8971 -1.0342 -1.9085\n" \
#    "H  5.2398  0.1213 -2.0263\n" \
#    "H  6.0816 -0.6698  0.1694\n" \
#    "H  4.8832 -1.7659  0.2457\n" \
#    "H -3.3202 -2.2775  0.3645\n" \
#    "H -5.6499  0.8241  0.0453"


# Ala-Gly-Ser-Glu   Total Energy = -1321.9687202425043
#
# geometry = "O  1.0081   -0.5551   -1.8575\n" \
#   "O -0.6789   -2.7971    0.2603\n" \
#   "O  5.4836   -2.3478   -0.3523\n" \
#   "O  3.9878   -2.4407    1.3520\n" \
#   "O -2.3256    1.3042    1.0021\n" \
#   "O  5.6406    3.0898    1.3743\n" \
#   "O  5.5752    2.8961   -0.8853\n" \
#   "O -5.2797   -0.9780   -0.5769\n" \
#   "N  2.2521   -0.8295    0.0829\n" \
#   "N -1.2590   -0.0770   -0.5337\n" \
#   "N -4.5827    1.2363   -0.5984\n" \
#   "N -7.8747   -0.1438   -0.1569\n" \
#   "C  3.5371   -0.9876   -0.5499\n" \
#   "C  4.2561    0.3485   -0.7723\n" \
#   "C -0.1287   -0.4409    0.2812\n" \
#   "C  1.0786   -0.6237   -0.6322\n" \
#   "C  4.6069    1.1378    0.4955\n" \
#   "C -0.4101   -1.6889    1.1117\n" \
#   "C  4.3332   -1.9811    0.2711\n" \
#   "C -2.2725    0.7719   -0.1037\n" \
#   "C -3.2967    1.0163   -1.1911\n" \
#   "C  5.3036    2.4532    0.2219\n" \
#   "C -6.7829    0.7201    0.2850\n" \
#   "C -5.4878    0.2088   -0.3332\n" \
#   "C -6.6835    0.7081    1.8053\n" \
#   "H  3.3628   -1.4575   -1.5263\n" \
#   "H  3.6334    0.9834   -1.4166\n" \
#   "H  5.1819    0.1669   -1.3337\n" \
#   "H  0.0877    0.3946    0.9585\n" \
#   "H  2.2209   -0.8149    1.0983\n" \
#   "H  3.6924    1.3604    1.0562\n" \
#   "H  5.2703    0.5378    1.1282\n" \
#   "H  0.4259   -1.9546    1.7660\n" \
#   "H -1.2920   -1.5393    1.7434\n" \
#   "H -1.3150   -0.4740   -1.4675\n" \
#   "H -3.3528    0.1682   -1.8806\n" \
#   "H -3.0128    1.9147   -1.7469\n" \
#   "H -4.8269    2.1883   -0.3415\n" \
#   "H -6.9822    1.7407   -0.0611\n" \
#   "H  0.1408   -3.0100   -0.2169\n" \
#   "H -6.4815   -0.2995    2.1869\n" \
#   "H -7.6095    1.0692    2.2656\n" \
#   "H -5.8654    1.3534    2.1442\n" \
#   "H  5.9905   -3.0034    0.1726\n" \
#   "H -8.7515    0.1741    0.2546\n" \
#   "H -7.9889   -0.0613   -1.1665\n" \
#   "H  6.1084    3.9346    1.2017\n"
#
# try:
#     psi4.core.set_output_file('output.dat', False)
#     psi4.set_memory('500 MB')
#
#     psi4.geometry(geometry)
#
#     psi4.energy('scf/cc-pvdz')
#     total_energy = psi4.core.get_variable('SCF TOTAL ENERGY')
#     print(total_energy)
# except RuntimeError:
#     pass

# ============================================

""" Normalize Matrix """

# matrix = [[ 1.2492,  1.1165, -0.4047],
# [ 2.5516, -2.2446, -0.9942],
# [ 0.5915, -1.1035, -1.1285],
# [ 1.7707,  1.2076, -1.7786],
# [-1.4105,  1.1507,  0.1821],
# [-0.7085, -0.1136,  0.3937],
# [ 0.7470,  0.0903,  0.0308],
# [-1.3345, -1.2000, -0.4702],
# [ 2.6493,  0.0652, -1.5860],
# [ 1.8114, -1.1318, -1.2310],
# [-0.7666, -0.3737,  1.4558],
# [-0.8580, -2.1695, -0.2878],
# [-2.4023, -1.3127, -0.2521],
# [-1.2248, -0.9797, -1.5384],
# [-2.3916,  1.0420,  0.4376],
# [-1.4071,  1.3875, -0.8099],
# [ 3.3505,  0.2709, -0.7726],
# [ 3.2005, -0.1396, -2.5078],
# [ 1.0947,  1.0073, -2.5146],
# [ 1.9913, -3.0108, -0.7463]]
#
# for i in range(len(matrix)):
#     matrix[i] = normalize_vector(matrix[i])
#
# for v in matrix:
#     print(v)

# ============================================

"""" Pearson Correlation """

# x_axis_reference = [1.2492,
# -1.4105,
# -0.7085,
# -1.3345,
#  0.7470,
# -0.7666,
# -0.8580,
# -2.4023,
# -1.2248,
# -2.3916,
# -1.4071,
#  2.5516,
#  0.5915,
#  1.7707,
#  2.6493,
#  1.8114,
#  3.3505,
#  3.2005,
#  1.0947,
#  1.9913]
#
#
#
# x_axis_result = [ 0.7678,
#  3.2928,
#  2.1762,
#  2.2573,
#  0.8780,
#  2.2275,
#  1.4277,
#  3.1926,
#  2.1960,
#  4.1724,
#  3.3017,
# -1.8033,
# -3.5624,
# -0.1435,
# -1.4706,
# -2.3923,
# -1.8824,
# -1.4058,
#  0.0278,
# -2.4226]
#
# results = np.corrcoef(x_axis_reference, x_axis_result)[0,1]
# print(results)


print("THE END")
sys.exit(0)
