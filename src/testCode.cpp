# include <vector>
# include <iostream>
# include <numeric>
# include <execution>
# include "genericSolverEA/Population.h"
# include "problemSpecificBackend/VrpRepXmlReader.h"


void swapVector(std::vector<int> &chromosome) {
    std::swap(chromosome[1], chromosome[10]);
}

void testVectorScope(std::vector<std::vector<int>> &population) {
    int idcs[4] = {1, 2, 3, 4};

    std::vector<std::vector<int>> offspring(4, std::vector<int>(3));

    int k = 1;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            offspring[i][j] = k++;
        }
    }

    for (int i = 0; i < 4; ++i) {
        //chromosomes[idcs[i]] = std::move(offspring[i]);
        //chromosomes[idcs[i]].assign(offspring[i].begin(), offspring[i].end());
        population[idcs[i]] = offspring[i];
    }
}

void testSumVector() {
    std::vector<int> lbs(10);
    for (int i = 0; i < lbs.size(); ++i) {
        lbs[i] = i;
    }
    int sumLBS = 0;
    for (auto n: lbs) {
        //sumLBS += n;
        int a = 1;
    }

    std::cout << sumLBS << std::endl;

    //SAMPLE
    std::vector<int> vec = {2, 4, 6, 8, 10, 12, 14, 16, 18};

//WITHOUT EXECUTION POLICY
    int sum = std::reduce(vec.begin(), vec.end());

//TAKING THE ADVANTAGE OF EXECUTION POLICIES
    int sum2 = std::reduce(std::execution::par, vec.begin(), vec.end());

    std::cout << "Without execution policy  " << sum << std::endl;
    std::cout << "With execution policy  " << sum2 << std::endl;

}

void testIntroToStdVector() {
    std::vector<int> vec(11);
    for (int i = 0; i < 11; ++i) {
        vec[i] = i;
    }
    swapVector(vec);
    for (int i = 0; i < 11; ++i) {
        std::cout << vec[i] + ", ";
    }
    std::cout << std::endl;

    std::vector<std::vector<int>> population(10, std::vector<int>(5));
    testVectorScope(population);

    std::vector<int> arrVectors[10];
}

void testPopulationInitPopulation() {

}

void testRandom() {
//    std::random_device rd; // obtain a random number from hardware
//    std::mt19937 gen(rd()); // seed the generator
    std::mt19937 gen{std::random_device{}()};
    std::uniform_int_distribution<> distr(25, 63); // define the range
    for (int n = 0; n < 40; ++n)
        std::cout << distr(gen) << ' '; // generate numbers
}

void test(std::vector<bool> &tr) {

    for (int i = 0; i < 10; ++i) {
        if (tr[i]) {
            std::cout << i << std::endl;
        }

    }
}

void testPopulationInit() {
    int populationSize = 10;
    int offspringSize = populationSize;
    int numberOfNodes = 10;

    std::vector<int> LBs(numberOfNodes);
    std::vector<int> UBs(numberOfNodes);
    std::fill(LBs.begin(), LBs.end(), 1);
    std::fill(UBs.begin(), UBs.end(), 1);
    LBs[0] = 2;
    UBs[0] = numberOfNodes + 1;

    Population population(populationSize, numberOfNodes, LBs, UBs);

    Population offspring(offspringSize);

}

class SimplePrinter {
public:
    void printVector(std::vector<int> &vec);
};

void SimplePrinter::printVector(std::vector<int> &vec) {
    for (auto &i: vec) {
        std::cout << i << std::endl;
    }
}

inline float euclideanDistance2D(std::pair<float, float> x, std::pair<float, float> y) {
    return std::sqrt(std::pow(x.first - y.first, 2) + std::pow((x.second - y.second), 2));
}

void testEuclideanDistance() {
    std::pair<float, float> x{40.0, 50.0};
    std::pair<float, float> y{45.0, 70.0};
    std::cout << euclideanDistance2D(x, y) << std::endl;
}

void testCrossover() {
//    VRPTWInstance vrptwInstance = createVRPTWInstance("../data/solomon-1987-c1/C101_025.xml");
//
//
//    GenericEASolver gs(&vrptwInstance, params);
//    GenericEASolver(&vrptwInstance, params).evolutionaryAlgorithm()
//
//
//    std::vector<int> chromosome1 = {0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 7, 0, 8, 0, 9, 0};
//    std::vector<int> chromosome2 = {0, 3, 2, 4, 0, 8, 5, 6, 1, 7, 9, 0};
//    printVec(chromosome1);
//    printVec(chromosome2);
//    std::vector<int> child = gs.crossover(chromosome2, chromosome1);
//    printVec(child);
}

void printVec(std::vector<int> &vec) {
    for (auto &i: vec) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

void testFitness() {
//    VRPTWInstance vrptwInstance = createVRPTWInstance("../data/solomon-1987-c1/C101_025.xml");
//    std::vector<int> perm = {0, 10, 11, 0, 22, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 2, 0,
//                             1, 0, 12, 0, 13, 0, 14, 0, 15, 0, 16, 0, 17, 0, 18, 0, 19, 0, 20,
//                             0, 21, 3, 23, 0, 24, 0, 25, 0};
//    std::cout << vrptwInstance.computeFitness(perm) << std::endl;
}

int main() {

    int depotID = 0;
    std::vector<int> goxChild = {0, 0,1, 2, 3, 0, 1, 0, 0, 1, 4, 0, 0};
    for (int i = 0; i < goxChild.size() - 1; ++i) {
        if (goxChild[i] == depotID && goxChild[i + 1] == depotID) {
            goxChild.erase(goxChild.begin() + i);
        }
    }

    std::vector<int> vec = {2, 3, 0, 1, 5};
    auto it = std::min(vec.begin(), vec.end());
    std::cout << *it << std::endl;

    std::vector<int> v;
    v.push_back(1);
    v.reserve(10);
    v.emplace_back(1);
    v.emplace_back(2);
    v.emplace_back(3);
    v.emplace_back(4);
    v.emplace_back(5);
    int sz = v.size();


    VrpRepXmlReader xmlReader("../data/solomon-1987-c1/C101_025.xml");

    std::vector<std::pair<float, float>> coords = xmlReader.getNodesCoordinates(10);

    int fs = xmlReader.getFleetSize();

    std::vector<float> st = xmlReader.getServiceTimes(10);
    int debugStop = 1;
}


// https://stackoverflow.com/questions/3221812/how-to-sum-up-elements-of-a-c-vector



/*
 *  OK
    std::random_device rd;
    std::mt19937 gen(rd());
    const std::uniform_int_distribution<> distr(10, 230);
    int randint = distr(gen);

    Has no members error
    std::mt19937 gen(std::random_device());
    const std::uniform_int_distribution<> distr(10, 230);
    int randint = distr(gen);


 A: https://stackoverflow.com/questions/19036141/vary-range-of-uniform-int-distribution
 */