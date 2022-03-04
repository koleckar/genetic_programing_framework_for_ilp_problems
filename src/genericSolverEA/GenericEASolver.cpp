# include "GenericEASolver.h"
# include <iostream>
# include <algorithm>
# include <random>
# include <cstdlib>
# include <cmath>
# include <numeric>
# include <limits>
# include <set>

using std::vector;


GenericEASolver::GenericEASolver(const Instance *instance,
                                 const Hyperparameters params) : instance(instance),
                                                                 params(params),
                                                                 population(params.populationSize,
                                                                            instance->numberOfNodes,
                                                                            instance->LBs,
                                                                            instance->UBs),
                                                                 offspring(params.offspringSize) {}


//////////// Helper Functions ////////////

inline bool equals(const char *str1, const char *str2) {
    return strcmp(str1, str2) == 0;
}

// todo: Q: in general, can there be the two same nodes next to each other in the permutation?
//  A: yes.
inline bool depotFreeNeighborhood(vector<int> &chromosome, int idx, int depotID) {
    if (chromosome[idx] == depotID) {
        return false;
    }
    if (idx > 0 && idx < chromosome.size() - 1) {
        return chromosome[idx + 1] != depotID && chromosome[idx - 1] != depotID;
    } else if (idx == 0) {
        return chromosome[idx + 1] != depotID;
    } else if (idx == chromosome.size() - 1) {
        return chromosome[idx - 1] != depotID;
    }
}

inline bool depotInsertionPossible(vector<int> &chromosome, int idx, int depotID) {
    if (idx == 0 || idx == chromosome.size() - 1) {
        return false;
    }
    if (chromosome[idx] == depotID) {
        return false;
    }
    return chromosome[idx - 1] != depotID;
}

inline bool genesCanBeSwapped(vector<int> &chromosome, int idx1, int idx2, int depotID = 0) {

    if (chromosome[idx1] != depotID && chromosome[idx2] != depotID) {
        return true;
    } else if (chromosome[idx1] == depotID && chromosome[idx2] != depotID) {
        if (idx1 == 0 || idx1 == chromosome.size() - 1) {
            return false;
        } else {
            return depotFreeNeighborhood(chromosome, idx2, depotID);
        }
    } else if (chromosome[idx1] != depotID && chromosome[idx2] == depotID) {
        if (idx2 == 0 || idx2 == chromosome.size() - 1) {
            return false;
        } else {
            return depotFreeNeighborhood(chromosome, idx1, depotID);
        }
    } else {
        return true;  // chromosome[idx1] == depotID && chromosome[idx2] == depotID
    }
}

inline int GenericEASolver::drawInt(int lb, int ub) {
    return std::uniform_int_distribution<int>{lb, ub}(randomNumberGenerator);
}

inline float GenericEASolver::drawFloat(float lb, float ub) {
    return std::uniform_real_distribution<float>{lb, ub}(randomNumberGenerator);
}

inline vector<int> kRandomIndicesFromRange(int first, int last, int k, std::mt19937 &randomNumberGenerator) {
    vector<int> rangeVector(last - first + 1); // both fist and last included.;
    std::iota(rangeVector.begin(), rangeVector.end(), first);

    std::vector<int> kRandomFromRange;
    kRandomFromRange.reserve(k);
    std::sample(rangeVector.begin(), rangeVector.end(),
                std::back_inserter(kRandomFromRange), k, randomNumberGenerator);
    return kRandomFromRange;
}

inline void swapIfFirstBigger(int &first, int &second) noexcept {
    if (first > second) {
        int tmp = first;
        first = second;
        second = tmp;
    }
}

inline void printStats(Population &population, int generation) {
    float bsf = std::numeric_limits<float>::max();
    int idx = 0;
    int idxBsf = 0;
    float f = -1;
    int feasibleCount = 0;
    int infeasibleCount = 0;
    float fSum = 0;
    int chromosomesLentghSum = 0;
    for (int i = 0; i < population.numberOfChromosomes; ++i) {
        chromosomesLentghSum += population.chromosomes[i].size();
        if (population.feasible[i]) {
            f = population.fitness[i];
            fSum += f;
            if (f < bsf) {
                bsf = f;
                idxBsf = idx;
            }
            feasibleCount++;
        } else {
            infeasibleCount++;
        }
        idx++;
    }
    std::cerr << generation << ": f* = " << bsf << "   cvs* = "
              << population.constraintsViolationSum[idxBsf]
              << "    fsbl* :" << population.feasible[idxBsf] << "   f_avg =" << fSum / feasibleCount << "   "
              << feasibleCount << "/"
              << population.numberOfChromosomes << "    avg_chr_len="
              << chromosomesLentghSum / population.numberOfChromosomes
              << std::endl;
}


//////////// EA Pipeline ////////////

void GenericEASolver::evolutionaryAlgorithm() {
    int generation = 0;

    initializePopulation(params.populationInitializationMode);

    calculateFitnessAndConstrainsViolationsSum(&population);

    bool iterate = true;
    while (iterate) {

        if (params.verbose) {
            if (generation % params.verbosePeriod == 0) {
                printStats(population, generation);
            }
        }

        calculateFitnessAndConstrainsViolationsSum(&population);

        breed();

        // todo: can you get along without repairing the offspring ?
        //repairOffspring();
        calculateFitnessAndConstrainsViolationsSum(&offspring);

        populationReplacement(params.replacementStrategyMode);

        generation++;
        if (generation++ >= params.generations) {
            iterate = false;
            std::cout << "EA stopping search." << std::endl;
        }
    }
}

void GenericEASolver::initializePopulation(const char *mode) {
    std::cout << "initializing population ..." << std::endl;

    if (equals(mode, "random-valid")) {
        vector<int> chromosome;
        bool feasible;
        int tries;
        for (int i = 0; i < population.numberOfChromosomes; ++i) {
            feasible = false;
            tries = 0;
            while (feasible == false && tries < 10) {
                chromosome = population.initializeChromosome(randomNumberGenerator);
                feasible = instance->isFeasible(chromosome);
                tries++;
            }
            population.chromosomes[i] = chromosome;
        }
    }
}

void GenericEASolver::breed() {

    vector<int> parentIndices = parentSelection(population.fitness,
                                                population.constraintsViolationSum, params.parentsSize);
    //vector<int> parentIndices = parentSelection(population.fitness, params.parentsSize);

    std::shuffle(parentIndices.begin(), parentIndices.end(), randomNumberGenerator);

    for (int i = 0; i < params.offspringSize; ++i) {
        // todo: wrong indexing ?
        int idxParent1 = parentIndices[drawInt(0, params.parentsSize - 1)];
        int idxParent2 = parentIndices[drawInt(0, params.parentsSize - 1)];
        offspring.chromosomes[i] = crossover(population.chromosomes[idxParent1],
                                             population.chromosomes[idxParent2]);
    }

    for (int i = 0; i < params.offspringSize; ++i) {
        if (drawFloat(0.0, 1.0) < params.mutationProb) {
            mutate(offspring.chromosomes[i], params.mutationMode);
        }
    }
    int deb = 1;
}

vector<int> GenericEASolver::stochasticRanking(vector<float> fitness, vector<float> scv, float pF = 0.45f) {
    /**
    "Stochastic Ranking" procedure
    for sorting of individuals in the population based on the objective function and constraint violation simultaneously
    as described in:
    https://www.researchgate.net/publication/3418601_Yao_X_Stochastic_ranking_for_constrained_evolutionary_optimization_IEEE_Trans_Evol_Comput_4_284-294

    Input: population, their evaluation = [f(x), sum(Gi(x))_i=<1, m+p> ] = [fitness, sum of constrains violations]
    Output: sorted population from best to worst one

    :param scv: summed constrained violations
    :param p_f: Probability of comparing two individuals only by their fitness. For both feasible ones p_f = 1.

    Note:
    -----------------------
     Use it with rank-based tournament selection within EA with generational replacement strategy.
    ----------------------
    */


    std::vector<int> indices;
    indices.reserve(fitness.size());
    for (int i = 0; i < fitness.size(); ++i) {
        indices.emplace_back(i);
    }

    for (int _ = 0; _ < fitness.size(); ++_) {

        bool noSwap = true;
        for (int j = 0; j < fitness.size() - 1; ++j) {

            float u = drawFloat(0, 1);
            if ((scv[indices[j]] == 0 && scv[indices[j + 1]] == 0) || u < pF) {
                if (fitness[indices[j]] > fitness[indices[j + 1]]) {
                    std::swap(indices[j], indices[j + 1]);
                    noSwap = false;
                }
            } else {
                if (scv[indices[j]] > scv[indices[j + 1]]) {
                    std::swap(indices[j], indices[j + 1]);
                    noSwap = false;
                }
            }
        }
        if (noSwap) {
            break;
        }
    }
    return indices;
}

vector<int> GenericEASolver::parentSelection(vector<float> &populationFitness,
                                             const int parentsSize, const char *mode) {
    std::vector<int> parentIndices;
    parentIndices.reserve(parentsSize);

    if (strcmp(mode, "tournament") == 0) {

        vector<float> decreasingProbabilities(params.tournamentSize);

        float prob = 0.5f;
        for (int i = 0; i < params.tournamentSize; ++i) {
            decreasingProbabilities[i] = prob * std::pow(1 - prob, i);
        }

        float sumForNormalization = 0;
        for (const auto &i: decreasingProbabilities) {
            sumForNormalization += i;
        }

        std::for_each(decreasingProbabilities.begin(), decreasingProbabilities.end(),
                      [sumForNormalization](float &prob) { prob /= sumForNormalization; });


        std::discrete_distribution<int> weightedDistribution
                (decreasingProbabilities.begin(), decreasingProbabilities.end());

        vector<int> tournamentIndices;
        for (int i = 0; i < parentsSize; ++i) {
            tournamentIndices = kRandomIndicesFromRange
                    (0, params.populationSize - 1, params.tournamentSize, randomNumberGenerator);

            std::sort(tournamentIndices.begin(), tournamentIndices.end(),
                      [&](int i, int j) { return population.fitness[i] > population.fitness[j]; });

            parentIndices.push_back(tournamentIndices[weightedDistribution(randomNumberGenerator)]);
        }
        return parentIndices;

    } else if (strcmp(mode, "by-fitness") == 0) {
        throw std::runtime_error("not implemented");
    } else {
        throw std::invalid_argument("wrong parent selection mode");
    }
}

vector<int> GenericEASolver::parentSelection(vector<float> &populationsFitness,
                                             vector<float> &sumOfConstraintsViolation,
                                             const int parentsSize) {
    std::vector<int> parentIndices;
    parentIndices.reserve(parentsSize);

    vector<int> populationRanks = stochasticRanking(populationsFitness, sumOfConstraintsViolation);
    // ranksAtIndices: key=specimen, value=rank, needed for the tournaments. The smaller the rank, the better.
    vector<int> ranksAtIndices(populationRanks.size());
    for (int rank = 0; rank < populationRanks.size(); ++rank) {
        int specimenIndex = populationRanks[rank];
        ranksAtIndices[specimenIndex] = rank;
    }

    vector<float> decreasingProbabilities(params.tournamentSize);
    float prob = 0.5f;
    for (int i = 0; i < params.tournamentSize; ++i) {
        decreasingProbabilities[i] = prob * std::pow(1 - prob, i);
    }
    float sumForNormalization = 0;
    for (const auto &i: decreasingProbabilities) {
        sumForNormalization += i;
    }
    std::for_each(decreasingProbabilities.begin(), decreasingProbabilities.end(),
                  [sumForNormalization](float &prob) { prob /= sumForNormalization; });

    std::discrete_distribution<int> weightedDistribution
            (decreasingProbabilities.begin(), decreasingProbabilities.end());

    vector<int> tournamentIndices;
    for (int i = 0; i < parentsSize; ++i) {
        tournamentIndices = kRandomIndicesFromRange
                (0, params.populationSize - 1, params.tournamentSize, randomNumberGenerator);

        // The smaller the rank, the better.
        std::sort(tournamentIndices.begin(), tournamentIndices.end(),
                  [&](int i, int j) { return ranksAtIndices[i] < ranksAtIndices[j]; });

        parentIndices.emplace_back(tournamentIndices[weightedDistribution(randomNumberGenerator)]);
    }
    return parentIndices;
}

vector<int> GenericEASolver::crossoverGOX(vector<int> &chromosome1, vector<int> &chromosome2) {
    int chromosomeLength = chromosome1.size(); // sizes do equal

    std::vector<int> child(chromosomeLength);

    vector<int> ordinals1(chromosomeLength); // at index i is ordinal of that gene at this position.
    vector<int> ordinals2(chromosomeLength);
    vector<int> geneCounts1(chromosomeLength);
    vector<int> geneCounts2(chromosomeLength);
    for (int i = 0; i < chromosome1.size(); ++i) {
        int gene1 = chromosome1[i];
        int gene2 = chromosome2[i];
        ordinals1[i] = ++geneCounts1[gene1];
        ordinals2[i] = ++geneCounts2[gene2];
    }

    int leftmostInsertedIdx = drawInt(0, chromosomeLength - 1);
    int rightmostInsertedIdx = drawInt(0, chromosomeLength - 1);
    swapIfFirstBigger(leftmostInsertedIdx, rightmostInsertedIdx);

    int insertionLength = rightmostInsertedIdx - leftmostInsertedIdx + 1;

    std::set<std::pair<int, int>> genesToRemove;
    for (int i = leftmostInsertedIdx; i <= rightmostInsertedIdx; ++i) {
        int gene = chromosome1[i];
        genesToRemove.emplace(gene, ordinals1[i]);
    }
    std::pair<int, int> firstInsertedGene
            (chromosome1[leftmostInsertedIdx], ordinals1[leftmostInsertedIdx]);

    // It is not known beforehand where the insertion in child will be.
    // First the prefix from chromosome 2 must be determined.
    // Then insert the genes from chromosome 1.
    // Then solve the suffix.

    // Solve prefix.
    int idx = 0;
    int prefixLength = 0;
    int chr2SuffixStartIndex;
    for (int i = 0; i < chromosomeLength; ++i) {
        int gene = chromosome2[i];
        int ordinal = ordinals2[i];
        if (firstInsertedGene.first == gene && firstInsertedGene.second == ordinal) {
            prefixLength = idx;
            chr2SuffixStartIndex = i;
            break;
        }
        if (genesToRemove.find({gene, ordinal}) != genesToRemove.end()) {
            continue;
        }
        child[idx++] = gene;
    }

    // Insert cut out part from chromosome 1.
    int j = leftmostInsertedIdx;
    for (int i = prefixLength; i < prefixLength + insertionLength; ++i) {
        child[i] = chromosome1[j++];
    }

    // Solve suffix.
    int childSuffixIdx = prefixLength + insertionLength;
    for (int i = chr2SuffixStartIndex; i < chromosomeLength; ++i) {
        int gene = chromosome2[i];
        int ordinal = ordinals2[i];

        if (genesToRemove.find({gene, ordinal}) != genesToRemove.end()) {
            continue;
        }
        child[childSuffixIdx] = gene;
        childSuffixIdx++;
    }

    /// todo--- Sanity check, debugging. -----
    vector<int> geneCount(population.numberOfGenes);
    for (int i = 0; i < child.size(); ++i) {
        int gene = child[i];
        geneCount[gene]++;
    }
    for (int i = 1; i < geneCount.size(); ++i) {
        if (geneCount[i] != 1) {
            throw std::runtime_error("wtf");
        }
    }
    ///-------------------------------------

    return child;
}

vector<int> GenericEASolver::crossover(vector<int> &chromosome1, vector<int> &chromosome2) {
    if (chromosome1.size() == chromosome2.size()) {

        vector<int> goxChild = crossoverGOX(chromosome1, chromosome2);


        // Repair child.
        for (int i = 0; i < goxChild.size() - 1; ++i) {
            if (goxChild[i] == depotID && goxChild[i + 1] == depotID) {
                goxChild.erase(goxChild.begin() + i);
                i--;
            }
        }
        if (goxChild[0] != depotID) {
            goxChild.insert(goxChild.begin(), depotID);
        }
        if (goxChild[goxChild.size() - 1] != depotID) {
            goxChild.push_back(depotID);
        }

        // todo: re-insert k random?
        // Reinsert depots, to match parent sizes.
        if (goxChild.size() < chromosome1.size()) {
            int diff = chromosome1.size() - goxChild.size();
            int c = 0;
            int idx = 0;
            int offset = 0;
            while (c != diff) {
                if (depotInsertionPossible(goxChild, idx, depotID)) {
                    goxChild.insert(goxChild.begin() + idx, depotID);
                    c++;
                }
                idx++;
            }
        }

        return goxChild;

    } else {

        vector<int> *longerChromosome = chromosome1.size() > chromosome2.size() ? &chromosome1 : &chromosome2;
        vector<int> *shorterChromosome = chromosome1.size() < chromosome2.size() ? &chromosome1 : &chromosome2;

        /// Match chromosome lengths.

        // Choose k depot stops at random. ( k = the difference in depot stops)
        int depotStops = longerChromosome->size() - instance->numberOfNodes + 1;
        int differenceInDepotStops = std::abs(int(chromosome2.size() - chromosome1.size()));
        int first = 2;  // cannot choose the first and last depot stop
        int last = depotStops - 1;
        int k = differenceInDepotStops;
        vector<int> removedDepotOrdinals = kRandomIndicesFromRange(first, last, k, randomNumberGenerator);

        // Remove the chosen depot stops from the longer chromosome.
        vector<int> shortenedLongerChromosome(shorterChromosome->size());
        int lastDepotOrdinal = 0;
        int idx = 0;
        vector<int> removedDepotIndices;
        removedDepotIndices.reserve(removedDepotOrdinals.size());
        for (int i = 0; i < longerChromosome->size(); ++i) {
            bool ok = true;
            if ((*longerChromosome)[i] == 0) {
                lastDepotOrdinal++;
                for (auto &depotOrdinalToRemove: removedDepotOrdinals) {
                    if (lastDepotOrdinal == depotOrdinalToRemove) {
                        ok = false;
                        removedDepotIndices.emplace_back(i);
                        break;
                    }
                }
            }
            if (ok == false) {
                continue;
            }
            shortenedLongerChromosome[idx++] = (*longerChromosome)[i];
        }

        vector<int> goxChild = crossoverGOX(shortenedLongerChromosome, *shorterChromosome);

        /// Repair child. (GOX can produce two same nodes next to each other)
        for (int i = 0; i < goxChild.size() - 1; ++i) {
            if (goxChild[i] == depotID && goxChild[i + 1] == depotID) {
                goxChild.erase(goxChild.begin() + i);
                i--;
            }
        }
        if (goxChild[goxChild.size() - 1] != depotID) {
            goxChild.push_back(depotID);
        }
        if (goxChild[0] != depotID) {
            goxChild.insert(goxChild.begin(), depotID);
        }

        /// Match length to shorter parent.
        if (goxChild.size() < shorterChromosome->size()) {
            int diff = shorterChromosome->size() - goxChild.size();
            int c = 0;
            int idx = 0;
            while (c != diff) {
                if (depotInsertionPossible(goxChild, idx, depotID)) {
                    goxChild.insert(goxChild.begin() + idx, depotID);
                    c++;
                }
                idx++;
            }
        }

        /// Repair at length to longer parent, with some lvl of randomness, metric is the block lengths.

        // Get sizes of blocks before and after removal,
        //  for each index get the change in blocks (tours) its removal makes.
        vector<std::pair<int, int>> blocksInLongerChromosomeMergedByDepotRemoval;
        blocksInLongerChromosomeMergedByDepotRemoval.reserve(differenceInDepotStops);
        for (int i = 0; i < differenceInDepotStops; ++i) {
            int depotIdx = removedDepotIndices[i];

            int leftIdx = depotIdx - 1;
            while ((*longerChromosome)[leftIdx] != depotID) {
                leftIdx--;
            }
            int leftBlockLength = depotIdx - leftIdx - 1;

            int rightIdx = depotIdx + 1;
            while ((*longerChromosome)[rightIdx] != depotID) {
                rightIdx++;
            }

            int rightBlockLength = rightIdx - depotIdx - 1;
            blocksInLongerChromosomeMergedByDepotRemoval.emplace_back(leftBlockLength, rightBlockLength);
        }


        // Get sizes of blocks in goxChild.
        vector<int> blockLengthsChild;
        vector<std::pair<int, int>> blockBordersChild;
        int maxPossibleBlocks = shorterChromosome->size() - instance->numberOfNodes;  // = Number of depot stops - 1
        blockLengthsChild.reserve(maxPossibleBlocks);
        blockBordersChild.reserve(maxPossibleBlocks);
        int depotIdx = 0;
        idx = 1;
        int blockLength = 0;
        while (depotIdx < goxChild.size() - 1) {
            if (goxChild[idx] == depotID) {
                blockLengthsChild.emplace_back(blockLength);
                int blockLeftIdx = idx - blockLength;
                int blockRightIdx = idx - 1;
                blockBordersChild.emplace_back(blockLeftIdx, blockRightIdx);
                depotIdx = idx;
                idx = depotIdx + 1;
                blockLength = 0;
            } else {
                idx++;
                blockLength++;
            }
        }

        // Choose at random number of block lengths to repair.
        int numberOfBlocksToRepair = drawInt(1, differenceInDepotStops);
        //int numberOfBlocksToRepair = differenceInDepotStops;

        // Choose blocks from goxChild that could be repaired.
        //  if more choices, choose at random
        //  if none of block lengths in child matches the merged blocks from longer chromosome terminate.
        int offset = 0; // Number of depots re-inserted so far.
        for (int i = 0; i < numberOfBlocksToRepair; ++i) {
            bool blockRepaired = false;
            for (int j = 0; j < blockLengthsChild.size(); ++j) {
                if (blockRepaired) {
                    break;
                }
                int childBlockLength = blockLengthsChild[j];
                for (int k = 0; k < blocksInLongerChromosomeMergedByDepotRemoval.size(); ++k) {
                    int mergedBlocksLength = blocksInLongerChromosomeMergedByDepotRemoval[k].first +
                                             blocksInLongerChromosomeMergedByDepotRemoval[k].second;
                    if (blockLength == mergedBlocksLength) {
                        int oneOfTheBlocksLengthBeforeMerge = blocksInLongerChromosomeMergedByDepotRemoval[k].second;
                        int idxOfRightBorderOfChildBlock = blockBordersChild[j].second;
                        int insertionIdx = idxOfRightBorderOfChildBlock + offset - oneOfTheBlocksLengthBeforeMerge;

                        goxChild.insert(goxChild.begin() + insertionIdx, depotID);

                        blockLengthsChild.erase(blockLengthsChild.begin() + j);
                        blockBordersChild.erase(blockBordersChild.begin() + j);
                        blocksInLongerChromosomeMergedByDepotRemoval.erase(
                                blocksInLongerChromosomeMergedByDepotRemoval.begin() + k);
                        blockRepaired = true;
                        offset += 1;
                        break;
                    }
                }
            }
            if (blockRepaired == false) { // If no repair was done, terminate.
                break;
            }
        }


        //todo --- Sanity check, debug ----
        for (int i = 0; i < goxChild.size() - 1; ++i) {
            if (goxChild[i] == depotID && goxChild[i + 1] == depotID) {
                throw std::runtime_error("wrong depot reinsertion");
            }
        }
        vector<int> geneCount(population.numberOfGenes);
        for (int i = 0; i < goxChild.size(); ++i) {
            int gene = goxChild[i];
            geneCount[gene]++;
        }
        for (int i = 1; i < geneCount.size(); ++i) {
            if (geneCount[i] != 1) {
                throw std::runtime_error("wtf");
            }
        }
        if (goxChild[0] != depotID || goxChild[goxChild.size() - 1] != depotID) {
            throw std::runtime_error("wtf");
        }
        //-------------------------------------

        return goxChild;
    }
}

void GenericEASolver::mutate(vector<int> &chromosome, const char *mode) {
    float u = drawFloat(0.1f, 1.0f);
    if (u < 0.3f) {
        mode = "swap2";
    } else if (0.3f <= u && u <= 0.5f) {
        mode = "delete-non-essential";
    } else {
        mode = "add-non-essential";
    }
    if (chromosome.size() < 31) {
        mode = "add-non-essential";
    }

    if (equals(mode, "swap2")) {
        int idx1 = drawInt(0, chromosome.size() - 1);
        int idx2 = drawInt(0, chromosome.size() - 1);
        if (genesCanBeSwapped(chromosome, idx1, idx2)) {
            std::swap(chromosome[idx1], chromosome[idx2]);
        }

    } else if (equals(mode, "delete-non-essential")) {
        int numberOfExtraDepotStops = chromosome.size() - instance->numberOfNodes - 1;
        int randomDepot;
        if (numberOfExtraDepotStops > 0) {
            randomDepot = drawInt(1, numberOfExtraDepotStops);
        } else {
            return;
        }
        int depotCounter = 0;
        for (int i = 2; i < chromosome.size() - 2; ++i) {
            if (chromosome[i] == depotID) {
                if (++depotCounter == randomDepot) {
                    chromosome.erase(chromosome.begin() + i);
                    return;
                }
            }
        }

    } else if (equals(mode, "add-non-essential")) {
        if (chromosome.size() == 2 * instance->numberOfNodes) {
            return;
        }
        // todo: make it more random, do not insert to first available every time.
        for (int i = 2; i < chromosome.size() - 1; ++i) {
            if (depotInsertionPossible(chromosome, i, depotID)) {
                chromosome.insert(chromosome.begin() + i, depotID);
                return;
            }
        }

    } else {
        throw std::invalid_argument("invalid mutation mode");
    }
}

void GenericEASolver::populationReplacement(const char *mode) {

    if (equals(mode, "generational")) {
        if (params.offspringSize == params.populationSize) {
            population.chromosomes = offspring.chromosomes;
        } else if (params.offspringSize > params.populationSize) {
            vector<int> indices(params.offspringSize);
            std::for_each(indices.begin(), indices.end(), [i = 0](int &x)mutable { x = i++; });

            vector<int> offspringRanks = stochasticRanking(offspring.fitness, offspring.constraintsViolationSum);
            for (int i = 0; i < params.populationSize; ++i) {
                population.chromosomes[i] = offspring.chromosomes[offspringRanks[i]];
            }
        } else {
            throw std::invalid_argument("offspring size < population size");
        }

    } else if (equals(mode, "steady-state")) {
        vector<vector<int>> chromosomes;
        vector<float> fitness;
        chromosomes = population.chromosomes;
        chromosomes.insert(chromosomes.end(), offspring.chromosomes.begin(), offspring.chromosomes.end());
        fitness = population.fitness;
        fitness.insert(fitness.end(), offspring.fitness.begin(), offspring.fitness.end());

        vector<vector<int>> survivors;
        survivors.reserve(params.populationSize);
        vector<int> alreadyChosen(params.populationSize + params.offspringSize, false);
        float min;
        float minIdx;
        for (int p = 0; p < params.populationSize; ++p) {
            min = std::numeric_limits<float>::max();
            minIdx = -1;

            for (int i = 0; i < chromosomes.size(); ++i) {
                if (alreadyChosen[i]) {
                    continue;
                }
                if (fitness[i] < min) {
                    min = fitness[i];
                    minIdx = i;
                }
            }
            alreadyChosen[minIdx] = true;
            survivors.emplace_back(chromosomes[minIdx]);

            if (survivors.size() == params.populationSize) {
                break;
            }
        }

        population.chromosomes = survivors;

    } else {
        throw std::invalid_argument("Invalid population replacement mode.");
    }
}

void GenericEASolver::calculateFitness(Population *specimen) {
    for (int i = 0; i < specimen->numberOfChromosomes; ++i) {
        specimen->fitness[i] = instance->computeFitness(specimen->chromosomes[i]);
        specimen->feasible[i] = instance->isFeasible(specimen->chromosomes[i]);
    }
}

void GenericEASolver::calculateFitnessAndConstrainsViolationsSum(Population *specimen, const char *mode) {
    for (int i = 0; i < specimen->numberOfChromosomes; ++i) {
        auto[chromosomeFitness, chromosomeCVS] =
        instance->computeFitnessAndConstraintsViolationSum(specimen->chromosomes[i]);

        specimen->fitness[i] = chromosomeFitness;
        specimen->constraintsViolationSum[i] = chromosomeCVS;
        float eps = 0.00001f;
        if (chromosomeCVS - eps > 0) {
            specimen->feasible[i] = false;
        } else {
            specimen->feasible[i] = true;
        }

        if (equals(mode, "death")) {
            specimen->fitness[i] += specimen->constraintsViolationSum[i] * 100000;
        }
    }
}



//vector<int> GenericEASolver::crossover_OLD_(vector<int> &chromosome1, vector<int> &chromosome2) {
//    if (chromosome1.size() == chromosome2.size()) {
//        return crossoverGOX(chromosome1, chromosome2);
//
//    } else {
//        vector<int> *longerChromosome = chromosome1.size() > chromosome2.size() ? &chromosome1 : &chromosome2;
//        vector<int> *shorterChromosome = chromosome1.size() < chromosome2.size() ? &chromosome1 : &chromosome2;
//
//        // todo:
//        //  in general GOX to function, there must be same count of every node in both chromosomes.
//        //  so in full generality, this condition must be repaired for all genes not just the depot gene.
//
//        // Choose at random the extra depots but from all depot stops.
//        // Choose k depot stops at random. k = the difference in depot stops.
//        int depotStops = longerChromosome->size() - instance->numberOfNodes + 1;
//        int differenceInDepotStops = std::abs(int(chromosome2.size() - chromosome1.size()));
//        int first = 2;  // cannot choose the first and last depot stop
//        int last = depotStops - 1;
//        int k = differenceInDepotStops;
//        vector<int> depotOrdinalsToRemove = kRandomIndicesFromRange(first, last, k, randomNumberGenerator);
//
//        // remove the chosen depot stops from the longer chromosome
//        vector<int> shortenedLongerChromosome(shorterChromosome->size());
//        int lastDepotOrdinal = 0;
//        int idx = 0;
//        vector<int> depotIndices;
//        depotIndices.reserve(depotOrdinalsToRemove.size());
//        for (int i = 0; i < longerChromosome->size(); ++i) {
//            bool ok = true;
//            if ((*longerChromosome)[i] == 0) {
//                lastDepotOrdinal++;
//                for (auto &depotOrdinalToRemove: depotOrdinalsToRemove) {
//                    if (lastDepotOrdinal == depotOrdinalToRemove) {
//                        ok = false;
//                        depotIndices.emplace_back(i);
//                        break;
//                    }
//                }
//            }
//            if (ok == false) {
//                continue;
//            }
//            shortenedLongerChromosome[idx++] = (*longerChromosome)[i];
//        }
//
//        vector<int> goxChild = crossoverGOX(shortenedLongerChromosome, *shorterChromosome);
//
//        // todo: Depot indices could be far over the indices of shorter chromosome.
//        //  lot of indices could be discarded !
//        depotIndices.erase(
//                std::remove_if(
//                        depotIndices.begin(),
//                        depotIndices.end(),
//                        [&](int i) {
//                            if (i > 0 && i < goxChild.size()) {
//                                return goxChild[i] == 0 || goxChild[i - 1] == 0;
//                            } else {
//                                return true;
//                            }
//                        }
//                ), depotIndices.end());
//
//        // todo: Do without resizing the vector k times.
//        // return zeros to goxChild;
//        int offset = 0;
//        for (auto &i: depotIndices) {
//            goxChild.insert(goxChild.begin() + i + offset++, 0);
//        }
//
//        //todo:
//        //  jak dava nejvic smysl doplnovat nuly?
//        //  asi co nejvic zachovat bloky.
//        //  predtim jsem mel extraDepotStops + 1  bloku takze bych mel zachovat aby
//        //  se delky po vraceni nul rovnaly jako pred odstranenim.
//
//        return goxChild;
//    }
//}


//void GenericEASolver::repairOffspring() {
//    for (auto &chromosome: offspring.chromosomes) {
//        repairChromosome(chromosome);
//    }
//}
//
//void GenericEASolver::repairChromosome(vector<int> &chromosome) {
//    // todo: use std algorithms ?
//    //      https://stackoverflow.com/questions/65528380/get-all-blocks-of-contiguous-same-numbers-for-an-array
//    vector<int> geneCounter(instance->numberOfNodes);
//    vector<int> repairedChromosome;
//
//    int startIdx, endIdx;
//    for (int i = 0; i < chromosome.size(); ++i) {
//        int gene = chromosome[i];
//        geneCounter[gene]++;
//
//        // skip the "invalid" block (consecutive depot stops)
//        if (gene == depotID) {
//            startIdx = i;
//            endIdx = i;
//            while (true) {
//                if (endIdx == chromosome.size() - 1) {
//                    break;
//                }
//                if (chromosome[endIdx + 1] == depotID) {
//                    endIdx++;
//                } else {
//                    break;
//                }
//            }
//            i = startIdx == endIdx ? startIdx : endIdx;
//        }
//        repairedChromosome.push_back(gene);
//    }
//    // todo: check correct gene counts.
//    if (repairedChromosome.size() < instance->numberOfNodes + 1) {
//        throw std::runtime_error("Some gene(s) are missing ?!?");
//    }
//    if (repairedChromosome[0] != depotID) {
//        repairedChromosome.insert(repairedChromosome.begin(), depotID);
//    }
//    if (repairedChromosome[repairedChromosome.size() - 1] != depotID) {
//        repairedChromosome.push_back(depotID);
//    }
//    chromosome = repairedChromosome;
//}
//

/// -------------------------------------------------------------------------------------
//inline vector<int> findGenesIndices(vector<int> &chromosome, int numberOfGenes) {
//    std::vector<int> indices(numberOfGenes);
//    for (int i = 0; i < chromosome.size(); ++i) {
//        indices[chromosome[i]] = i;
//    }
//    return indices;
//}
//
//inline void swapIfFirstBigger(int &i, int &j) {
//    if (i > j) {
//        int tmp = j;
//        j = i;
//        i = tmp;
//    }
//}

//vector<int> GenericEASolver::crossover(vector<int> &chromosome) {
//    vector<int> child(chromosome.size());
//
//    int numberOfExtraDepotStops = chromosome.size() - instance->numberOfNodes - 1;
//    // invarianty: pozice of block of genes mezi depoty v chromozomu
//    // ne invarianty: dva po sobe jdouci geny, fitness(g1g2) != (g2g1) - time(g1->g2) == time(g2->g1)
//    // tak casova okna mohou byt ruzna.
//    if (numberOfExtraDepotStops == 0) {
//
//    } else if (numberOfExtraDepotStops == 1) {
//
//    }
//    // swap two regions
//    // the regions between zeros.
//    return child;
//}

//vector<int> GenericEASolver::crossoverExperiment(vector<int> &chromosome1, vector<int> &chromosome2) {
//    // 1. find block of sames its max length;
//    // vyhledavani vzoru v textu ! najdi nejdelsi schodu podretezce k
//    vector<int> chromosome2GeneIndices = findGenesIndices(chromosome2, instance->numberOfNodes);
//
//    vector<int> blocksOfSameGenesOffsets(instance->numberOfNodes); // [startGene] = endGene (offset)
//
//    int bestID = -1;
//    int bestLength = -1;
//    for (int i = 0; i < chromosome1.size(); ++i) {
//        if (chromosome1[i] == 0) {
//            continue;
//        }
//        int idx1 = i;
//        int idx2 = chromosome2GeneIndices[i];
//        bool search = true;
//        while (search) {
//            if (idx1 == chromosome1.size() || idx2 == chromosome2.size()) {
//                search = false;
//            }
//            if (chromosome1[idx1] == 0) {
//                idx1 += 1;
//            }
//            if (chromosome2[idx2] == 0) {
//                idx2 += 1;
//            }
//
//            if (chromosome1[idx1] == chromosome2[idx2]) {
//                idx1++;
//                idx2++;
//            } else {
//                search = false;
//            }
//        }
//        int blockLength = idx1 - i;
//        if (blockLength > bestLength) {
//            bestID = i;
//            bestLength = blockLength;
//        }
//        blocksOfSameGenesOffsets[i] = idx1;
//    }
//    // swap the blocks
//    return vector<int>(10);
//}
//
//vector<int> GenericEASolver::crossoverWithModes(vector<int> &chromosome1, vector<int> &chromosome2, const char *mode) {
//    // todo: if the chromosomes lengths equals use gox.
//
//    if (equals(mode, "GOX")) {
//        return crossoverGOX(chromosome1, chromosome2);
//
//    } else if (equals(mode, "exG")) {
//        //todo:
//    } else if (equals(mode, "SVLC")) {
//        //todo:
//    } else if (equals(mode, "GOX-experimental")) {
//        // 1. equal lenghts by removing / adding zeros
//        // 2. call GOX
//        // 3. add zeros back to the position
//
//    } else if (equals(mode, "experiment")) {
//        crossoverExperiment(chromosome1, chromosome2);
//    } else {
//        throw std::invalid_argument("invalid crossoverWithModes mode");
//    }
//}