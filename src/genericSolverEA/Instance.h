# ifndef SEMESTRAL_WORK_INSTANCE_H
# define SEMESTRAL_WORK_INSTANCE_H

# include <vector>
# include <memory>
# include <string>

/**
 Interface (abstract class) representing one instance of some optimization problem (eg. CVRP, VRPTW, TSP,...)
 representable by permutation with repetition and variable length.
 */
class Instance {
protected:
    Instance::Instance(int numberOfNodes, std::vector<int> LBs, std::vector<int> UBs, std::string instanceName)
            : numberOfNodes(numberOfNodes), LBs(LBs), UBs(UBs), instanceName(instanceName) {}

public:
    const int numberOfNodes;
    const std::vector<int> LBs, UBs;
    const std::string instanceName;

    virtual float computeFitness(const std::vector<int> &permutation) const = 0;

    virtual std::pair<float, float> computeFitnessAndConstraintsViolationSum(const std::vector<int> &permutation) const = 0;

    virtual bool isFeasible(const std::vector<int> &permutation) const = 0;
};


#endif //SEMESTRAL_WORK_INSTANCE_H
