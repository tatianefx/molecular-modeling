from modules.chemistry.amino_acid import AminoAcid
from modules.common.helper import rotation_euler
from modules.common.helper import calculates_xyz_to_rotate
from functools import reduce
import psi4


class Individual:

    def __init__(self, chromosome: [AminoAcid]):
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

        # print(geometry + '\n')
        total_energy = 0
        try:
            psi4.core.set_output_file('output.dat', False)
            psi4.set_memory('500 MB')

            psi4.geometry(geometry)

            psi4.energy('scf/cc-pvdz')
            total_energy = psi4.core.get_variable('SCF TOTAL ENERGY')
        except RuntimeError:
            # Iterations did not converge.
            pass

        return total_energy

    def mutate(self):
        for i in range(len(self.chromosome)):
            first_part = self.chromosome[:i+1]
            second_part = self.chromosome[i+1:]

            items = reduce(lambda x, y: x+y, second_part)
            part_to_rotate = []
            for item in items:
                for atom in item.atoms:
                    part_to_rotate = atom.position

            # calculates resulting vector
            v = calculates_xyz_to_rotate(first_part[-1].atoms[-1].position, part_to_rotate[0])

            # move positions
            new_part = rotation_euler(part_to_rotate, v)

            # update values
            for j in range(len(second_part)):
                for k in range(len(second_part[j])):
                    (second_part[j])[k] = new_part[0]
                    del new_part[0]

            # update origin
            self.chromosome[i + 1:] = second_part

        # calculate fitness
        self.fitness = self.__calculates_fitness()
        print('Fitness: ' + str(self.fitness))

    def __write3d(self):
        str = ""
        for gene in self.chromosome:
            # TODO transform chromosome to 3D string
            return str
        return str
