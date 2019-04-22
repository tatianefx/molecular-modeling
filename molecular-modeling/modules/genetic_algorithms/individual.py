from modules.genetic_algorithms.gene import Gene
from modules.common.helper import rotation_euler
import numpy as np


class Individual:

    def __init__(self, chromosome: [Gene]):
        self.chromosome = chromosome
        self.fitness = self.__calculates_fitness()

    def __calculates_fitness(self):
        # TODO Calculates the Gibbs free energy.

        return 0.0

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