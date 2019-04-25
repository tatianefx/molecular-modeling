from modules.genetic_algorithms.gene import Gene
from modules.common.helper import rotation_euler
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
        geometry = []

        peptide = self.chromosome
        for item in peptide:
            for atom in item.atoms:
                geometry.append(atom.position)

        # split at random position

        first_part = geometry[0:int(len(geometry)/2)]
        second_part = geometry[int(len(geometry)/2):len(geometry)]

        # move positions
        new_part = rotation_euler(second_part, [0, 0, 1])

        # update geometry
        geometry.clear()
        geometry += first_part
        geometry += new_part

        for item in peptide:
            for i in range(0, len(item.atoms)):
                item.atoms[i].position = geometry[i]

        # calculate fitness
        self.fitness = self.__calculates_fitness()
        print("Fitness: " + str(self.fitness))

    def __write3d(self):
        str = ""
        for gene in self.chromosome:
            # TODO transform chromosome to 3D string
            return str
        return str