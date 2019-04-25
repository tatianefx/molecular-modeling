#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "tatianefx"

from modules.genetic_algorithms.individual import Individual
from modules.common.config import POPULATION_SIZE, GENERATIONS
import random


class GeneticAlgorithm:

    the_best_individual: Individual

    def __init__(self, individual: Individual):
        self.first_individual = individual
        self.population = self.__generate_population()
        self.the_best_individual = self.__run_genetic_algorithm()

    def __run_genetic_algorithm(self):
        for i in range(GENERATIONS):
            self.__select_population()
            self.__crossover()
            self.__mutate()
        self.__select_population()
        return self.population[0]

    def __generate_population(self):
        population = []
        for i in range(POPULATION_SIZE*2):
            individual = Individual(self.first_individual)
            individual.mutate()
            population.append(individual)
        return population

    def __select_population(self):
        self.population.sort(key=lambda x: x.fitness, reverse=True)
        del self.population[POPULATION_SIZE:]

    def __crossover(self):
        offsprings = []
        for i in range(0, POPULATION_SIZE, 2):
            offspring = self.__uniform_crossover(self.population[i], self.population[i+1])
            offsprings.append(offspring)
        self.population = self.population + offsprings

    def __mutate(self):
        for i in range(len(self.population)):
            self.population[i] = self.population[i].mutate()

    @staticmethod
    def __uniform_crossover(parent_a: Individual, parent_b: Individual):
        genes = []
        for i in range(len(parent_a.chromosome)):
            choice = random.randrange(2)
            if choice == 0:
                genes.append(parent_a.chromosome[i])
            else:
                genes.append(parent_b.chromosome[i])
        return Individual(genes)

