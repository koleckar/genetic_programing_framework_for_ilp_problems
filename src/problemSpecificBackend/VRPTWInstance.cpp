#include "VRPTWInstance.h"

#include <iostream>
#include <cmath>

using std::vector;
using std::pair;


//TODO:
//  * Add variants of hyper-parameters to fitness calculation.
//  * merge the functions doing almost the same thing.


VRPTWInstance::VRPTWInstance(int numberOfNodes, vector<int> LBs, vector<int> UBs, std::string instanceName) :
        Instance(numberOfNodes, LBs, UBs, instanceName),
        timeMatrix(numberOfNodes, vector<float>(numberOfNodes)) {}


inline float euclideanDistance2D(std::pair<float, float> &x, pair<float, float> &y) {
    return std::sqrt(std::pow(x.first - y.first, 2) + std::pow((x.second - y.second), 2));
}

vector<vector<float>> VRPTWInstance::createDistanceMatrixFromCoordinates(vector<pair<float, float>> &coordinates) {

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

void VRPTWInstance::createTimeMatrixFromDistanceMatrix(vector<vector<float>> &distanceMatrix) {
    float normalizationFactor = 1;

    for (int i = 0; i < numberOfNodes; ++i) {
        for (int j = 0; j < numberOfNodes; ++j) {
            timeMatrix[i][j] = distanceMatrix[i][j] * normalizationFactor;
        }
    }
}

void VRPTWInstance::initializeDataStructures(VrpRepXmlReader &xmlReader) {
    vector<pair<float, float>> coordinates = xmlReader.getNodesCoordinates(numberOfNodes);
    vector<vector<float>> distanceMatrix = createDistanceMatrixFromCoordinates(coordinates);
    createTimeMatrixFromDistanceMatrix(distanceMatrix);

    serviceTimes = xmlReader.getServiceTimes(numberOfNodes);

    timeWindows = xmlReader.getTimeWindows(numberOfNodes);

    demandsNodes = xmlReader.getNodeDemands(numberOfNodes);

    capacitiesVehicles = std::vector<float>(xmlReader.getFleetSize(), xmlReader.getVehicleCapacity());
}


inline constexpr float deathPenalty() noexcept {
    return 1000000;
}

pair<float, float> VRPTWInstance::computeOneTour(const vector<int> &permutation, int start, int end,
                                                 const char *mode) const {
    float time = 0;
    float vehicleCapacity = capacitiesVehicles[1];
    float penalization = 0;

    int fromID, toID;
    for (int i = start; i < end; ++i) {
        fromID = permutation[i];
        toID = permutation[i + 1];
        time += timeMatrix[fromID][toID];

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

        if (time < timeWindows[toID].first) {
            time = timeWindows[toID].first;
        } else if (time > timeWindows[toID].second) {
            if (strcmp(mode, "death") == 0) {
                penalization += deathPenalty();
            }
            if (strcmp(mode, "cvs") == 0) {
                penalization += time - timeWindows[toID].second;
            }

        } else {
            // constraints ok.
        }

        time += serviceTimes[toID];
    }
    return pair<float, float>{time, penalization};
}

float VRPTWInstance::computeFitness(const std::vector<int> &permutation) const {
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


pair<float, float> VRPTWInstance::computeFitnessAndConstraintsViolationSum(const std::vector<int> &permutation) const {
    /**
     Computes fitness but no penalties are added for the constraints violation.
     Instead all constraints violations are summed with no penalties.
        Eg. time=115, tw=<60, 80>  -> constraintsViolationSum += 35
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

// todo: switch for homogeneous/heterogeneous fleet ? or consider only homogeneous.
float VRPTWInstance::isTourFeasible(const vector<int> &permutation, int start, int end) const {
    bool isFeasible = true;
    int vehicleCapacity = capacitiesVehicles[1];

    float time = 0;
    float penalization = 0;

    int fromID, toID;
    for (int i = start; i < end; ++i) {
        fromID = permutation[i];
        toID = permutation[i + 1];
        time += timeMatrix[fromID][toID];

        if (toID == depotID) {
            break;
        }

        vehicleCapacity -= demandsNodes[toID];
        if (vehicleCapacity < 0) {
            return false;
        }

        if (time < timeWindows[toID].first) {
            time = timeWindows[toID].first;
        } else if (time > timeWindows[toID].second) {
            return false;
        } else {
            // constraints ok.
        }

        time += serviceTimes[toID];
    }

    return true;
}

// todo: can be done with std algorithms? (increase readability)
bool VRPTWInstance::isFeasible(const vector<int> &permutation) const {
    int length = permutation.size();
    std::vector<bool> visitsSatisfied(numberOfNodes);
    std::vector<int> visits(numberOfNodes);

    int node;
    for (int i = 0; i < length; ++i) {
        node = permutation[i];
        visits[node]++;
    }
    for (int i = 1; i < numberOfNodes; ++i) {
        if (visits[i] != 1) {
            return false;
        }
    }

    if (permutation[0] != depotID || permutation[length - 1] != depotID) {
        return false;
    }

    for (int i = 1; i < length; ++i) {
        if (permutation[i] == depotID) {
            if (permutation[i - 1] == depotID) {
                return false;
            }
        }
    }

    bool feasible;
    int depotIdx1 = 0, depotIdx2 = 0;
    for (int i = 0; i < permutation.size(); ++i) {
        if (permutation[i] == depotID) {
            depotIdx1 = depotIdx2;
            depotIdx2 = i;
            feasible = isTourFeasible(permutation, depotIdx1, depotIdx2);
            if (feasible == false) {
                return false;
            }
        }
    }
    return true;
}

