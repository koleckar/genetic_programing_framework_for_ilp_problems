# Genetic programing framework for ILP problem
<h4>Genetic programing framework for various ILP problems with instances representable as permutation with repetition</h4>
______________________________________________________________________________________________________________________________________________________
<h5>Work motivation</h5>
This work aims at creating metaheuristics optimization solver/framework (mof)[1] for optimization problems with solutions that can be meaningfully represented as permutations with repetition and variable length (pwravl). 
The permutations and permutations with repetition being subset of prwavl. 
Example: Given a set {0, 1, 2, 3, 4},  some of solutions using this representation could be:

[0, 1, 2, 3, 4],  [0, 1, 0, 2, 4]  and  [0, 0, 0, 1, 2, 3, 4, 1, 2, 3].  
(Nothing about validity of such solutions for any problem is said, they stand here solely as example of the pwravl representation.)


Solving of large or hard instances of NP-hard problems (i.e. their optimization version) is often intractable even by the best commercial MILP solvers (eg. Gurobi).
The MILP solvers (and ILP formulation) strength is the ease of implementation and also if the solver is able to solve the problem instance, it is solved to the global optimum.

Even though lots of optimization problems are subjected to intense research and advanced methods that solve hard instances of ‘standardized’ datasets -even to the optimum- can be found in literature, they are often non-trivial, very problem specific and hard to implement. Also the palette of metaheuristics usable for addressing hard problems is notably large. 
Single-solution (e.g. iterated local search, tabu-search) vs. population-based (genetic algorithms, ant-colony optimization), the combinations of population and single solution approaches in the form of memetic algorithms, etc.

Even for someone proficient in this domain however, it can be intricate to choose and implement the most suitable one.
This work aims at being somewhere in the middle. Given any problem representable by pwravl, the only thing that user has to implement is the fitness (cost) function for the given problem instance, following a simple predefined interface. The solver should be able to tackle even the hard problems. The generality and ease of use comes however with a cost. Such approach cannot compete with the problem-tailored heuristics and some gap for the optimum relative to optimum found by problem-tailored heuristics is presumable.

