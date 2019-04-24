from modules.genetic_algorithms.gene import Gene
from modules.common.helper import rotation_euler
import numpy as np
import psi4


class Individual:

    def __init__(self, chromosome: [Gene]):
        self.chromosome = chromosome
        self.fitness = self.__calculates_fitness()

    def __calculates_fitness(self):
        geometry = ""

        for item in self.chromosome:
            for atom in item.atoms:
                if geometry != "":
                    geometry = geometry + '\n'
                geometry = geometry + atom.acronym.value + ' '
                geometry = geometry + str("%.4f" % atom.position[0]) + ' '
                geometry = geometry + str("%.4f" % atom.position[1]) + ' '
                geometry = geometry + str("%.4f" % atom.position[2])

        psi4.core.set_output_file('output.dat', False)
        psi4.set_memory('500 MB')

        psi4.geometry(geometry)

        psi4.energy('scf/cc-pvdz')
        total_energy = psi4.core.get_variable('SCF TOTAL ENERGY')

        return total_energy

    def mutate(self):
        # get 3D position
        geometry = np.ndarray()

        peptide = self.chromosome
        for item in peptide:
            for atom in item.atoms:
                geometry += atom.position

        # split at random position
        parts = np.split(geometry, 2)

        # move positions
        new_part = rotation_euler(parts[1], [0, 0, 1])

        # update geometry
        geometry = parts[0] + new_part

        for item in peptide:
            for atom in item.atoms:
                atom.position = geometry[atom.index(atom.position)]

        # calculate fitness
        self.fitness = self.__calculates_fitness()

    def __write3d(self):
        str = ""
        for gene in self.chromosome:
            # TODO transform chromosome to 3D string
            return str
        return str