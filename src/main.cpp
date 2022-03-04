# include <iostream>
# include <numeric>

# include "genericSolverEA/Population.h"
# include "genericSolverEA/Solution.h"
# include "genericSolverEA/GenericEASolver.h"

# include "problemSpecificBackend/VRPTWInstance.h"
# include "problemSpecificBackend/CVRPInstance.h"
# include "problemSpecificBackend/VrpRepXmlReader.h"

# include "../lib/tinyxml/tinyxml.h"


// TODO : (main pipeline)
//  * (DONE) set hyper-parameters (check their validity)
//  * (DONE) run EA (print intermediate results)
//  * (DONE) choose best from result population -> create Solution instance
//  * (DONE) return the Solution - write it to xml file.

// todo:  (other tasks)
//  * be compliant with the integer rounding VRP datasets propose?
//  * adjust depots re-insertion in crossover()
//  * add time limit ?
//  * (DONE) initialize better the population, now most of chromosomes are infeasible.
//  * start with only feasible solutions
//      -> for larger datasets valid-random is computationally heavy, used some kind of informed init instead?
//          or the simples valid with full length permutations?
//  * (DONE) use ea for constrained optimization problem -> stochastic ranking
//  * move helper free functions to some tool source file
//  * compare the speed increase after using threads
//  * add more info to solution xml file
//      -> hyperparameters used, info about solution: waiting times, number of depot stops, etc...
//  * (DONE) save hyperparameters
//  * write time elapsed to stats
//  * depotID as member variable? depotID exists only if all lbs==ubs except one having ubs>lbs?
//
// todo:  (special crossover)
//  * return two children?
//  * shortening the chromosomes to fast?
//  * maybe not solved in full generality?
//      check: some for loops maybe written assuming that [0] [-1] is depotID but for eg. TSP there would be no depotID !!!
//      holds tho for permutations(wr) with constant length


// todo: (evaluation)
//  * compare to:
//      + ILP



Hyperparameters createHyperparameters() {
    int populationSize = 1000;
    int offspringSize = 5000;
    int parentsSize = (int) std::floor(populationSize * 0.8f);
    int generations = 20000;
    const char *populationInitializationMode = "random-valid";
    const char *replacementStrategyMode = "generational";
    const char *parentSelectionMode = "";
    int tournamentSize = (int) std::floor(float(populationSize) / 12.0);
    const char *crossoverMode = "";
    const char *mutationMode = "delete-non-essential"; // swap2, delete-non-essential, add-non-essential
    float mutationProb = 0.4f;
    bool verbose = true;
    int verbosePeriod = (int) generations / 200
            ;
    float stochasticRankingProb = 0.5f;

    return Hyperparameters{populationSize, offspringSize, parentsSize, generations, populationInitializationMode,
                           replacementStrategyMode, parentSelectionMode, tournamentSize, crossoverMode, mutationMode,
                           mutationProb, verbose, verbosePeriod, stochasticRankingProb};
}

bool areHyperparametersValid(const Hyperparameters *params) {
    if (params->offspringSize < params->populationSize) {
        //throw std::invalid_argument("offspringSize < populationSize");
        std::cerr << "offspringSize < populationSize" << std::endl;
        return false;
    }

    if (params->tournamentSize < 3) {
        //throw std::invalid_argument("tournament size must be >= 3.");
        std::cerr << "tournament size must be >= 3." << std::endl;
        return false;
    }

    if (params->parentsSize < (int) params->populationSize / 3) {
        //throw std::invalid_argument("parentSize must be at least third of populationSize");
        std::cerr << "parentSize must be at least third of populationSize" << std::endl;
        return false;
    }

    return true;
}


VRPTWInstance createVRPTWInstance(const char *filePath) {
    VrpRepXmlReader xmlReader(filePath);

    int numberOfNodes = xmlReader.getNumberOfNodes();

    std::vector<int> LBs(numberOfNodes);
    std::vector<int> UBs(numberOfNodes);
    std::fill(LBs.begin(), LBs.end(), 1);
    std::fill(UBs.begin(), UBs.end(), 1);
    LBs[0] = 2;
    UBs[0] = numberOfNodes + 1;

    VRPTWInstance vrptwInstance(numberOfNodes, LBs, UBs, xmlReader.datasetName);

    vrptwInstance.initializeDataStructures(xmlReader);

    return vrptwInstance;
}

CVRPInstance createCVRPInstance(const char *filePath) {
    VrpRepXmlReader xmlReader(filePath);

    int numberOfNodes = xmlReader.getNumberOfNodes();

    std::vector<int> LBs(numberOfNodes);
    std::vector<int> UBs(numberOfNodes);
    std::fill(LBs.begin(), LBs.end(), 1);
    std::fill(UBs.begin(), UBs.end(), 1);
    LBs[0] = 2;
    UBs[0] = numberOfNodes + 1;

    CVRPInstance cvrpInstance(numberOfNodes, LBs, UBs, xmlReader.datasetName);

    cvrpInstance.initializeDataStructures(xmlReader);

    return cvrpInstance;
}


Solution getOptimalSolutionFromPopulation(GenericEASolver &gs) {
    auto iteratorOpt = std::min_element(gs.population.fitness.begin(), gs.population.fitness.end());
    int indexOpt = iteratorOpt - gs.population.fitness.begin();

    float fitnessOpt = gs.population.fitness.at(indexOpt);
    std::vector<int> chromosomeOpt = gs.population.chromosomes.at(indexOpt);
    bool isFeasibleOpt = gs.population.feasible.at(indexOpt);
    float cvsOpt = gs.population.constraintsViolationSum.at(indexOpt);

    return Solution(chromosomeOpt, fitnessOpt, isFeasibleOpt, cvsOpt);
}

inline std::string stdVectorToString(std::vector<int> &vec) noexcept {
    std::string s = "";
    for (auto &i: vec) {
        s += std::to_string(i);
        s += " ";
    }
    s.pop_back();
    return s;
}

void saveSolutionToXML(Solution solution, Hyperparameters params, const char *filePath, const char *instanceName) {
    TiXmlDocument xmlDocument;

    TiXmlDeclaration *declaration = new TiXmlDeclaration("1.0", "", "");
    xmlDocument.LinkEndChild(declaration);

    TiXmlElement *root = new TiXmlElement("solution");
    root->SetAttribute("dataset_instance_name", instanceName);
    xmlDocument.LinkEndChild(root);

    TiXmlElement *solutionElem = new TiXmlElement("permutation");
    TiXmlElement *fitnessElem = new TiXmlElement("fitness");
    TiXmlElement *cvsElem = new TiXmlElement("constraints_violation_sum");
    TiXmlElement *feasibleElem = new TiXmlElement("feasible");

    std::string permutationStr = stdVectorToString(solution.permutation);
    solutionElem->LinkEndChild(new TiXmlText(permutationStr.c_str()));
    fitnessElem->LinkEndChild(new TiXmlText(std::to_string(solution.fitness).c_str()));
    cvsElem->LinkEndChild(new TiXmlText(std::to_string(solution.cvs).c_str()));
    feasibleElem->LinkEndChild(new TiXmlText(std::to_string(solution.isFeasible).c_str()));


    TiXmlElement *hyperparamsElem = new TiXmlElement("hyperparameters");

    TiXmlElement *populationSizeElem = new TiXmlElement("population_size");
    TiXmlElement *offspringSizeElem = new TiXmlElement("offspring_size");
    TiXmlElement *parentsSizeElem = new TiXmlElement("parents_size");
    TiXmlElement *generationsElem = new TiXmlElement("generations");
    TiXmlElement *populationInitializationModeElem = new TiXmlElement("population_initialization");
    TiXmlElement *parentSelectionModeElem = new TiXmlElement("parent_selection");
    TiXmlElement *tournamentSizeElem = new TiXmlElement("tournament_size");
    TiXmlElement *crossoverModeElem = new TiXmlElement("crossover");
    TiXmlElement *mutationModeElem = new TiXmlElement("mutation");
    TiXmlElement *mutationProbElem = new TiXmlElement("mutation_probability");
    TiXmlElement *stochasticRankingProbElem = new TiXmlElement("stochastic_ranking_probability");
    TiXmlElement *replacementStrategyModeElem = new TiXmlElement("replacement_strategy");

    populationSizeElem->LinkEndChild(new TiXmlText(std::to_string(params.populationSize).c_str()));
    offspringSizeElem->LinkEndChild(new TiXmlText(std::to_string(params.offspringSize).c_str()));
    parentsSizeElem->LinkEndChild(new TiXmlText(std::to_string(params.parentsSize).c_str()));
    generationsElem->LinkEndChild(new TiXmlText(std::to_string(params.generations).c_str()));
    populationInitializationModeElem->LinkEndChild(new TiXmlText(params.populationInitializationMode));
    parentSelectionModeElem->LinkEndChild(new TiXmlText(params.parentSelectionMode));
    tournamentSizeElem->LinkEndChild(new TiXmlText(params.crossoverMode));
    crossoverModeElem->LinkEndChild(new TiXmlText(params.crossoverMode));
    mutationModeElem->LinkEndChild(new TiXmlText(params.mutationMode));
    mutationProbElem->LinkEndChild(new TiXmlText(std::to_string(params.mutationProb).c_str()));
    stochasticRankingProbElem
            ->LinkEndChild(new TiXmlText(std::to_string(params.stochasticRankingProb).c_str()));
    replacementStrategyModeElem->LinkEndChild(new TiXmlText(params.replacementStrategyMode));

    hyperparamsElem->LinkEndChild(populationSizeElem);
    hyperparamsElem->LinkEndChild(offspringSizeElem);
    hyperparamsElem->LinkEndChild(parentsSizeElem);
    hyperparamsElem->LinkEndChild(generationsElem);
    hyperparamsElem->LinkEndChild(populationInitializationModeElem);
    hyperparamsElem->LinkEndChild(replacementStrategyModeElem);
    hyperparamsElem->LinkEndChild(tournamentSizeElem);
    hyperparamsElem->LinkEndChild(crossoverModeElem);
    hyperparamsElem->LinkEndChild(mutationProbElem);
    hyperparamsElem->LinkEndChild(stochasticRankingProbElem);

    root->LinkEndChild(solutionElem);
    root->LinkEndChild(fitnessElem);
    root->LinkEndChild(cvsElem);
    root->LinkEndChild(feasibleElem);
    root->LinkEndChild(hyperparamsElem);

    xmlDocument.SaveFile(filePath);
    std::cout << "solution saved" << " " << filePath << std::endl;
}


int main() {

    Hyperparameters params = createHyperparameters();
    if (areHyperparametersValid(&params) == false) {
        std::cerr << "Wrong hyperparameter(s) passed." << std::endl;
        return -1;
    }

    VRPTWInstance vrptwInstance = createVRPTWInstance("../data/solomon-1987-c1/C101_025.xml");
    GenericEASolver gs(&vrptwInstance, params);
    //CVRPInstance cvrpInstance = createCVRPInstance("../data/uchoa-et-al-2014/X-n1001-k43.xml");
    //GenericEASolver gs(&cvrpInstance, params);
    gs.evolutionaryAlgorithm();

//    Solution optimalSolution = getOptimalSolutionFromPopulation(gs);
//    std::string fileName = "../data-solutions/" + vrptwInstance.instanceName + ".xml";
//    saveSolutionToXML(optimalSolution, params, fileName.c_str(), vrptwInstance.instanceName.c_str());

//    Solution optimalSolution = getOptimalSolutionFromPopulation(gs);
//    std::string fileName = "../data-solutions/" + cvrpInstance.instanceName + ".xml";
//    saveSolutionToXML(optimalSolution, params, fileName.c_str(), cvrpInstance.instanceName.c_str());

    return 0;
}
