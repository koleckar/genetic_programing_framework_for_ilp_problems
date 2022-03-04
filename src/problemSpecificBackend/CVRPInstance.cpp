#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#include "CVRPInstance.h"

using std::vector;
using std::pair;


//TODO:
//  * Add variants of hyper-parameters to fitness calculation.
//  * merge the functions doing almost the same thing.


CVRPInstance::CVRPInstance(int numberOfNodes, vector<int> LBs, vector<int> UBs, std::string instanceName) :
        Instance(numberOfNodes, LBs, UBs, instanceName) {}


inline float euclideanDistance2D(std::pair<float, float> &x, pair<float, float> &y) {
    return std::sqrt(std::pow(x.first - y.first, 2) + std::pow((x.second - y.second), 2));
}

vector<vector<float>> CVRPInstance::createDistanceMatrixFromCoordinates(vector<pair<float, float>> &coordinates) {

    vector<vector<float>> distanceMatrix
            (coordinates.size(), vector<float>(coordinates.size()));

    int i = 0, j = 0;
    for (auto &coord1: coordinates) {
        for (auto &coord2: coordinates) {
            if (i == j) {
                distanceMatrix[i][j] = 0;
            } else {
                float distance = euclideanDistance2D(coord1, coord2);
                distanceMatrix[i][j] = distance;
                distanceMatrix[j][i] = distance;
            }
            j++;
        }
        i++;
        j = 0;
    }

    return distanceMatrix;
}

void CVRPInstance::initializeDataStructures(VrpRepXmlReader &xmlReader) {
    vector<pair<float, float>> coordinates = xmlReader.getNodesCoordinates(numberOfNodes);
    distanceMatrix = createDistanceMatrixFromCoordinates(coordinates);

    demandsNodes = xmlReader.getNodeDemands(numberOfNodes);

    // unlimited number of vehicles
    capacitiesVehicles = std::vector<float>(1, xmlReader.getVehicleCapacity());
}


inline constexpr float deathPenalty() noexcept {
    return 1000000;
}

pair<float, float> CVRPInstance::computeOneTour(const vector<int> &permutation, int start, int end,
                                                const char *mode) const {
    float distance = 0;
    float vehicleCapacity = capacitiesVehicles[0];
    float penalization = 0;

    int fromID, toID;
    for (int i = start; i < end; ++i) {
        fromID = permutation[i];
        toID = permutation[i + 1];
        distance += distanceMatrix[fromID][toID];

        if (toID == depotID) {
            break;
        }

        if (vehicleCapacity - demandsNodes[toID] < 0) {
            if (strcmp(mode, "death") == 0) {
                penalization += deathPenalty();
            }
            if (strcmp(mode, "cvs") == 0) {
                penalization += vehicleCapacity < 0 ? demandsNodes[toID] : demandsNodes[toID] - vehicleCapacity;
            }
        }
        vehicleCapacity -= demandsNodes[toID];
    }

    return pair<float, float>{distance, penalization};
}

float CVRPInstance::computeFitness(const std::vector<int> &permutation) const {
    /**
     For for constraints violation death penalty is added to the fitness. (simplest solution to constraints violation)
     */
    float fitness = 0;

    int depotIdx1 = 0, depotIdx2 = 0;
    for (int i = 0; i < permutation.size(); ++i) {
        if (permutation[i] == depotID) {
            depotIdx1 = depotIdx2;
            depotIdx2 = i;
            auto[tourFitness, tourPenalization] = computeOneTour(permutation, depotIdx1, depotIdx2);
            fitness += tourFitness + tourPenalization;
        }
    }

    return fitness;
}


pair<float, float> CVRPInstance::computeFitnessAndConstraintsViolationSum(const std::vector<int> &permutation) const {
    /**
     Computes fitness but no penalties are added for the constraints violation.
     Instead all constraints violations are summed with no penalties.
        Eg. capacity=15 demand=35 -> constraintsViolationSum += 20
     cvs used by "stochastic ranking" or as second objective in eg. NSGA2
     */
    float fitness = 0;
    float cvs = 0;

    int depotIdx1 = 0, depotIdx2 = 0;
    for (int i = 0; i < permutation.size(); ++i) {
        if (permutation[i] == depotID) {
            depotIdx1 = depotIdx2;
            depotIdx2 = i;
            auto[tourFitness, tourPenalization] = computeOneTour(permutation, depotIdx1, depotIdx2);
            fitness += tourFitness;
            cvs += tourPenalization;
        }
    }

    return pair<float, float>{fitness, cvs};
}

bool CVRPInstance::isFeasible(const std::vector<int> &permutation) const {
    return true;
}

