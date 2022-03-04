# include "Population.h"
# include <numeric>
# include <iostream>
# include <set>

Population::Population(int populationSize) : chromosomes(populationSize),
                                             fitness(populationSize),
                                             constraintsViolationSum(populationSize),
                                             feasible(populationSize),
                                             numberOfChromosomes(populationSize),
                                             numberOfGenes(-1) {
    // Constructor for offspring, only initializing member variables.
}

Population::Population(int populationSize,
                       int numberOfGenes,
                       std::vector<int> genesLBs,
                       std::vector<int> genesUBs,
                       const char *populationInitMode) : chromosomes(populationSize),
                                                         fitness(populationSize),
                                                         constraintsViolationSum(populationSize),
                                                         feasible(populationSize),
                                                         numberOfChromosomes(populationSize),
                                                         numberOfGenes(numberOfGenes),
                                                         genesLBs(genesLBs),
                                                         genesUBs(genesUBs) {
    // initializePopulation(populationInitMode);
}


inline int drawInt(int lb, int ub, std::mt19937 &randomNumberGenerator) {
    return std::uniform_int_distribution<int>{lb, ub}(randomNumberGenerator);
}

inline void solveOneGene(std::vector<int> &genes, std::vector<int> &chromosome, std::mt19937 &randomNumberGenerator) {
    int idx = drawInt(0, genes.size() - 1, randomNumberGenerator);
    int gene = genes.at(idx);
    chromosome.push_back(gene);
    genes.erase(genes.begin() + idx);
}

std::vector<int> Population::initializeChromosome(std::mt19937 &randomNumberGenerator) {
    /**
     Initialization for vrp problems. Depot is added first and last.
     At each cycle, a weighted coin is flipped, to add node or depot.
     Weight of adding a node decreases with a number of momentarily consecutive nodes.
     */
    std::vector<int> chromosome;
    chromosome.reserve(numberOfGenes);
    std::vector<int> genes;
    generate_n(inserter(genes, genes.begin()), numberOfGenes - 1, [i = 1]() mutable { return i++; });

    chromosome.push_back(0);
    solveOneGene(genes, chromosome, randomNumberGenerator);

    int consecutiveEssentials = 2;
    float prob = 0.5f;
    float probDepot, probNode;
    while (genes.empty() == false) {
        probNode = prob * std::pow(1 - prob, consecutiveEssentials + 1) ;
        probDepot = 1 - probNode;
        std::discrete_distribution<int> depotOrNode{probDepot, probNode};
        int weightedCoinFlip = depotOrNode(randomNumberGenerator);
        if (weightedCoinFlip == 0) {
            chromosome.push_back(0);
            solveOneGene(genes, chromosome, randomNumberGenerator);
            consecutiveEssentials = 1;

        } else {
            solveOneGene(genes, chromosome, randomNumberGenerator);
            consecutiveEssentials++;
        }
    }

    if (chromosome[chromosome.size() - 1] != 0) {
        chromosome.push_back(0);
    }

    return chromosome;
}


