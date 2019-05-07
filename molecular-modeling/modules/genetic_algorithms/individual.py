from modules.chemistry.amino_acid import AminoAcid
from modules.common.helper import rotation_euler
from modules.common.helper import calculates_xyz_to_rotate
from modules.common.helper import calculates_angle_a_b_c
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
            pass

        return total_energy

    def mutate(self):
        second_part = []

        for atom in self.chromosome[1].atoms:
            second_part.append(atom.position)

        a = self.chromosome[0].atoms[-2].position
        b = self.chromosome[0].atoms[-1].position
        c = second_part[0]
        before = calculates_angle_a_b_c(a, b, c)
        print('Angle before:' + str(before))

        # calculates resulting vector
        v = calculates_xyz_to_rotate(self.chromosome[0].atoms[-2].position, self.chromosome[0].atoms[-1].position)

        # move positions
        new_part = rotation_euler(second_part, v)

        a = self.chromosome[0].atoms[-2].position
        b = self.chromosome[0].atoms[-1].position
        c = new_part[0]
        after = calculates_angle_a_b_c(a, b, c)
        print('Angle after:' + str(after))

        for i in range(len(self.chromosome[1].atoms)):
            self.chromosome[1].atoms[i].position = new_part[i]

        # calculate fitness
        self.fitness = self.__calculates_fitness()
        print('Fitness: ' + str(self.fitness) + '\n')

    def __write3d(self):
        str = ""
        for gene in self.chromosome:
            # TODO transform chromosome to 3D string
            return str
        return str