#ifndef SEMESTRAL_WORK_CVRPINSTANCE_H
#define SEMESTRAL_WORK_CVRPINSTANCE_H

# include "VrpRepXmlReader.h"
# include "../genericSolverEA/Instance.h"

class CVRPInstance : public Instance {

    std::vector<float> demandsNodes;
    std::vector<float> capacitiesVehicles;
    std::vector<std::vector<float>> distanceMatrix;

    const int depotID = 0;  // Without loss of generality.

    std::vector<std::vector<float>> CVRPInstance::createDistanceMatrixFromCoordinates
            (std::vector<std::pair<float, float>> &coordinates);

    std::pair<float, float> CVRPInstance::computeOneTour(const std::vector<int> &permutation, int start, int end,
                                                         const char *mode = "cvs") const;

    float CVRPInstance::isTourFeasible(const std::vector<int> &permutation, int start, int end) const;

public:
    CVRPInstance(int numberOfNodes, std::vector<int> LBs, std::vector<int> UBs, std::string instanceName);

    float computeFitness(const std::vector<int> &permutation) const override;

    std::pair<float, float>
    computeFitnessAndConstraintsViolationSum(const std::vector<int> &permutation) const override;

    void initializeDataStructures(VrpRepXmlReader &xmlReader);

    bool isFeasible(const std::vector<int> &permutation) const override;

};


#endif //SEMESTRAL_WORK_CVRPINSTANCE_H
