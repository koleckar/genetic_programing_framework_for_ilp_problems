# include "ILPModel.h"
# include <gurobi_c++.h>
# include "../problemSpecificBackend/VrpRepXmlReader.h"
# include <vector>
# include <cmath>
# include <limits>

using std::vector;
using std::pair;


inline float euclideanDistance2D(std::pair<float, float> &x, pair<float, float> &y) {
    return std::sqrt(std::pow(x.first - y.first, 2) + std::pow((x.second - y.second), 2));
}


vector<vector<int>> createDistanceMatrixFromCoordinates(vector<pair<float, float>> &&coordinates) {

    vector<vector<int>> distanceMatrix
            (coordinates.size(), vector<int>(coordinates.size()));

    int i = 0, j = 0;
    for (auto &coord1: coordinates) {
        for (auto &coord2: coordinates) {
            if (i == j) {
                distanceMatrix[i][j] = 0;
            } else {
                int distance = round(euclideanDistance2D(coord1, coord2));

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

void twoIndexVehicleFlowFormulation(int N, int Q, int K, vector<vector<int>> distanceMatrix,
                                    vector<int> demands, vector<pair<int, int>> timeWindows) {
    /// ILP two-index flow formulation for CVRPTW, as in https://arxiv.org/pdf/1606.01935.pdf .
    GRBEnv env;
    GRBModel model(env);
    model.set(GRB_StringAttr_ModelName, "VRP-TW ILP model");

    vector<vector<GRBVar>> x;
    x.reserve((N + 1) * (N + 1));// binary, x_ij == 1 iff route goes from i->j
    for (int i = 0; i < N + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            std::string name = "x_" + std::to_string(i) + std::to_string(j);
            int objCoeffs;
            if (i == N) {
                objCoeffs = distanceMatrix[0][j];
            } else if (j == N) {
                objCoeffs = distanceMatrix[i][0];
            } else if (i == N && j == N) {
                objCoeffs = distanceMatrix[0][0];
            } else {
                objCoeffs = distanceMatrix[i][j];
            }
            x.emplace_back(model.addVar(0, 1, objCoeffs, GRB_BINARY, name));
        }
    }

    vector<GRBVar> y;
    y.reserve(N);
    for (int i = 0; i < N + 1; ++i) {
        std::string name = "y_i" + std::to_string(i);
        if (i == N) {
            y.emplace_back(model.addVar(demands[0], Q, 0, GRB_CONTINUOUS, name.c_str()));
        } else {
            y.emplace_back(model.addVar(demands[i], Q, 0, GRB_CONTINUOUS, name.c_str()));
        }
    }

    vector<GRBVar> w;
    w.reserve(N);
    for (int i = 0; i < N + 1; ++i) {
        std::string name = "w" + std::to_string(i);
        if (i == N) {
            y.emplace_back(model.addVar(
                    timeWindows[0].first, timeWindows[0].second, 0, GRB_CONTINUOUS, name.c_str()));
        } else {
            y.emplace_back(model.addVar(
                    timeWindows[i].first, timeWindows[i].second, 0, GRB_CONTINUOUS, name.c_str()));
        }
    }

    // Constraint: all customers are visited exactly once .
    for (int j = 1; j < N; ++j) {
        GRBLinExpr expr = 0;
        for (int i = 0; i < N; ++i) {
            if (i == j) {
                continue;
            }
            expr += x[i][j];
            std::string s = "customer " + std::to_string(i) + " visited exactly once";
            model.addConstr(expr == 1, s.c_str());
        }
    }

    // Correct flow of the vehicles.
    for (int i = 0; i < N; ++i) {

    }

    // Maximum number of vehicles leaving depot. == K means we allow self loop for the depot. eg- x[0][n] = 1
    GRBLinExpr expr = 0;
    for (int j = 0; j < N; ++j) {
        expr += x[0][j];
    }
    model.addConstr(expr <= K);

    // objective coefficients are set during the creation of the decision variables above.
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);


    // Solve the model.
    model.optimize();

    // Print the objective
    // and the values of the decision variables in the solution.
//    std::cout << "Optimal objective: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
//    std::cout << "x: " << x.get(GRB_DoubleAttr_X) << " ";
//    std::cout << "y: " << y.get(GRB_DoubleAttr_X) << std::endl;
//    std::cout << "w: " << y.get(GRB_DoubleAttr_X) << std::endl;
}

void threeIndexVehicleFlowFormulation(int N, int Q, int K, vector<vector<int>> distanceMatrix,
                                      vector<int> demands, vector<pair<int, int>> timeWindows) {
    /// as given in https://reader.elsevier.com/reader/sd/pii/S1018364710000297?token=9FADA6554ECCE12A5E12D35BD4A5B2ADD681E2213839F403847CB643836E778BF8671D023E24E2EF5CA3BB97BA423623&originRegion=eu-west-1&originCreation=20220114122755

}

// todo: will model of a LARGE instance fit to stack?
int main(int argc, char *argv[]) {

    // todo: code would be cleaner if first datastructures made N+1.

    /// Load data. All values are integers, if not they are rounded.
    VrpRepXmlReader vrpReader("../data/solomon-1987-c1/C101_025.xml");

    int N = vrpReader.getNumberOfNodes(); // Representation is: start at depot 0, return to depot N+1.

    int Q = vrpReader.getVehicleCapacity();

    int K = vrpReader.getFleetSize();

    vector<vector<int>> distanceMatrix  // Distances are rounded to the closest int.
            = createDistanceMatrixFromCoordinates(vrpReader.getNodesCoordinates(N));

    vector<int> demands;
    demands.reserve(N);
    for (auto f: vrpReader.getNodeDemands(N)) {
        demands.emplace_back((int) f);
    }

    vector<std::pair<int, int>> timeWindows;
    timeWindows.reserve(N);
    for (auto[f1, f2]: vrpReader.getTimeWindows(N)) {
        timeWindows.emplace_back((int) f1, (int) f2);
    }
    timeWindows[0] = std::pair<int, int>{0, std::numeric_limits<int>::max()};

    twoIndexVehicleFlowFormulation(N, Q, K, distanceMatrix, demands, timeWindows);

    return 0;
}


// ILP formulation cubic in variables, not used.
// arc set:
// (0, i);
// (i, j): 1) a_i + t_ij <= b_j , t_ij= time from i->j + serviceTime_i , (i, n+1)
//         2) d_i + d_j <= q
// (i, n+1);
// all values are integers.
// triangle inequality on cost holds.
// start_v_1 = a_v_1
// start_v_i = max {s_v_(i-1) + t[v_(i-1),v_i], a_v_i}