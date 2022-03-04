#ifndef SEMESTRAL_WORK_POPULATION_H
#define SEMESTRAL_WORK_POPULATION_H

#include <vector>
#include <random>


/**
 Class representing the population of solutions of the evolution algorithm (EA) solver .
 Genes could represent eg. the nodes in some solution of some optimization problem solved by the EA eg. TSP, CVRP, ...
    representable by the permutation with repetition and variable length (PWRVL).
 Chromosome is the PWRVL.
 */
class Population {

public:
    std::vector<std::vector<int>> chromosomes;
    std::vector<float> fitness;
    std::vector<float> constraintsViolationSum;
    std::vector<bool> feasible;

    const int numberOfGenes;
    const int numberOfChromosomes;
    const std::vector<int> genesLBs;
    const std::vector<int> genesUBs;


    // Constructor for offspring without initializing the chromosomes.
    Population(int populationSize);

    Population(int populationSize, int numberOfGenes, std::vector<int> genesLBs, std::vector<int> genesUBs,
               const char *populationInitMode = "valid-random");

    std::vector<int> initializeChromosome(std::mt19937 &randomNumberGenerator);

};

#endif //SEMESTRAL_WORK_POPULATION_H
