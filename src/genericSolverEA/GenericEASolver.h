#ifndef SEMESTRAL_WORK_GENERICEASOLVER_H
#define SEMESTRAL_WORK_GENERICEASOLVER_H

# include "Instance.h"
# include "Population.h"
# include <cstdlib>
# include <random>

struct Hyperparameters {
    int populationSize;
    int offspringSize;
    int parentsSize;
    int generations;
    const char *populationInitializationMode;
    const char *replacementStrategyMode;
    const char *parentSelectionMode;
    int tournamentSize;
    const char *crossoverMode;
    const char *mutationMode;
    float mutationProb;
    bool verbose;
    int verbosePeriod;
    float stochasticRankingProb;
    int timeout;
};

class GenericEASolver {
    const Instance *instance;
    const Hyperparameters params;
    int timeout;
    const int depotID = 0;

    std::mt19937 randomNumberGenerator{std::random_device{}()};

    void calculateFitness(Population *population);

    void calculateFitnessAndConstrainsViolationsSum(Population *population, const char *mode = "csv");

    void initializePopulation(const char *mode);

    void breed();

    std::vector<int> parentSelection(std::vector<float> &populationFitness,
                                     const int parentsSize, const char *mode = "tournament");

    std::vector<int> parentSelection(std::vector<float> &populationsFitness,
                                     std::vector<float> &sumOfConstraintsViolation,
                                     const int parentSize);

    std::vector<int> crossover(std::vector<int> &chromosome1, std::vector<int> &chromosome2);

    std::vector<int> crossoverGOX(std::vector<int> &chromosome1, std::vector<int> &chromosome2);

    void mutate(std::vector<int> &chromosome, const char *mode = "swap2");

    void populationReplacement(const char *mode = "generational");

    std::vector<int> stochasticRanking(std::vector<float> fitness, std::vector<float> scv, float pF);

    int GenericEASolver::drawInt(int lb, int ub);

    float GenericEASolver::drawFloat(float lb, float ub);

public:

    Population population;
    Population offspring;

    GenericEASolver(const Instance *instance, const Hyperparameters params);

    void evolutionaryAlgorithm();

};


#endif //SEMESTRAL_WORK_GENERICEASOLVER_H
