from modules.chemistry.amino_acid import AminoAcid
from modules.common.helper import rotation_euler
from modules.common.helper import calculates_xyz_to_rotate
from functools import reduce
import psi4
from random import randint


class Individual:

    def __init__(self, chromosome: [AminoAcid]):
        self.chromosome = chromosome
        self.fitness = self.__calculates_fitness()
        self.geometry = ''

    def __calculates_fitness(self):
        geometry = ""
        random_number = randint(1, 100)

        for item in self.chromosome:
            for atom in item.atoms:
                if geometry != "":
                    geometry = geometry + '\n'
                geometry = geometry + atom.acronym.value + ' '
                geometry = geometry + str("%.4f" % atom.position[0]) + ' '
                geometry = geometry + str("%.4f" % atom.position[1]) + ' '
                geometry = geometry + str("%.4f" % atom.position[2])

        # print(geometry + '\n')
        self.geometry = geometry
        total_energy = 0
        try:
            psi4.core.set_output_file('output'+str(random_number)+'.dat', False)
            psi4.set_memory('500 MB')

            psi4.geometry(geometry)

            psi4.energy('scf/cc-pvdz')
            total_energy = psi4.core.get_variable('SCF TOTAL ENERGY')
        except RuntimeError:
            # Iterations did not converge.
            pass

        return total_energy

    def mutate(self):
        if len(self.chromosome) == 1:
            return

        for i in range(1, len(self.chromosome)):
            first_part: list = self.chromosome[:i]
            second_part: list = self.chromosome[i:]

            part_to_rotate = []
            if len(second_part) == 1:
                items = reduce(lambda x, y: x + y, second_part)
                for atom in items.atoms:
                    part_to_rotate.append(atom.position)
            else:
                items = reduce(lambda x, y: x.atoms + y.atoms, second_part)
                for atom in items:
                    part_to_rotate.append(atom.position)

            # calculates resulting vector
            v = calculates_xyz_to_rotate(first_part[-1].atoms[-1].position, part_to_rotate[0])

            # move positions
            new_part = rotation_euler(part_to_rotate, v)

            # update values by reference
            for j in range(len(second_part)):
                for k in range(len(second_part[j].atoms)):
                    second_part[j].atoms[k].position = new_part[0]
                    del new_part[0]

        # calculate fitness
        self.fitness = self.__calculates_fitness()
        print('Fitness: ' + str(self.fitness))

    def __write3d(self):
        str = ""
        for gene in self.chromosome:
            # TODO transform chromosome to 3D string
            return str
        return str
