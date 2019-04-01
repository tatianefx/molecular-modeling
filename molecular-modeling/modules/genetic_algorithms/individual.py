from modules.genetic_algorithms.gene import Gene


class Individual:

    def __init__(self, chromosome: [Gene]):
        self.chromosome = chromosome
        self.fitness = self.__calculates_fitness()

    def __calculates_fitness(self):
        # TODO Calculates the Gibbs free energy.

        return 0.0

    def mutate(self):
        # TODO
        self.fitness = self.__calculates_fitness()

    def __write3d(self):
        str = ""
        for gene in self.chromosome:
            # TODO transform chromosome to 3D string
            return str
        return str