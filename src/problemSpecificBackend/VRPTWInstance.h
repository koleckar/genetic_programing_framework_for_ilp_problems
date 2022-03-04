#ifndef SEMESTRAL_WORK_VRPTW_H
#define SEMESTRAL_WORK_VRPTW_H

# include "VrpRepXmlReader.h"
# include "../genericSolverEA/Instance.h"


// todo: Read depotID from the dataset (even though it's always ==0)?

/**
Class representing instance of some VRP-TW optimization problem.

Problem hyper-parameters:
                * time matrix - time cost of traveling between the nodes.
                * time windows (hard or soft)
                * fleet: homogeneous, heterogeneous
                * vehicle: capacity (inf/bounded), max_travel_time
                * customer: demand, service_time (0, >0)

The objective is to minimize the time while satisfying all the instance constraints.
 */
class VRPTWInstance : public Instance {


    std::vector<std::pair<float, float>> timeWindows;
    std::vector<float> demandsNodes;
    std::vector<float> capacitiesVehicles;
    std::vector<float> serviceTimes;

    const int depotID = 0;  // Without loss of generality.

    void createTimeMatrixFromDistanceMatrix(std::vector<std::vector<float>> &distanceMatrix);

    std::vector<std::vector<float>> VRPTWInstance::createDistanceMatrixFromCoordinates
            (std::vector<std::pair<float, float>> &coordinates);

    std::pair<float, float> VRPTWInstance::computeOneTour(const std::vector<int> &permutation, int start, int end,
                                                          const char *mode = "cvs") const;

    float VRPTWInstance::isTourFeasible(const std::vector<int> &permutation, int start, int end) const;

public:
    std::vector<std::vector<float>> timeMatrix;

    VRPTWInstance(int numberOfNodes, std::vector<int> LBs, std::vector<int> UBs, std::string instanceName);

    float computeFitness(const std::vector<int> &permutation) const override;

    std::pair<float, float> computeFitnessAndConstraintsViolationSum(const std::vector<int> &permutation) const override;

    bool isFeasible(const std::vector<int> &permutation) const override;

    void initializeDataStructures(VrpRepXmlReader &xmlReader);

};

#endif //SEMESTRAL_WORK_VRPTW_H
