#ifndef SEMESTRAL_WORK_SOLUTION_H
#define SEMESTRAL_WORK_SOLUTION_H

# include <limits>
# include <vector>


// todo: is 'Solution' good name for Best Found Solution ?

/**
 Class representing the best found permutation solution to some optimization problem (eg. CVRPTW, TSP, Job-Shop...)
 representable by permutation with repetition and variable length.
 Best in means of minimal fitness value, if possible chosen in population from the feasible ones.
 */
class Solution {
public:
    int fitness;
    int cvs;
    bool isFeasible;
    std::vector<int> permutation;

    Solution(std::vector<int> &permutation, int fitness, bool isFeasible) :
            permutation(permutation), fitness(fitness), isFeasible(isFeasible) {};

    Solution(std::vector<int> &permutation, int fitness, int cvs, bool isFeasible) :
            permutation(permutation), fitness(fitness), cvs(cvs), isFeasible(isFeasible) {};

};

#endif //SEMESTRAL_WORK_SOLUTION_H
