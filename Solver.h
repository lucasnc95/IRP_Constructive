#ifndef SOLVER_H
#define SOLVER_H

#include "IRP.h"
#include "Solution.h"
#include "Vehicle.h"
//#include <Route.h>
#include <random>

class Solver {
public:
    const IRP& irp;
    std::mt19937 rng;

    Solver(const IRP& irp);

    Solution solve(int dmax);
    void updateInventory(int period, int customerId, std::vector<std::vector<int>>& currentInventory, int deliveryAmount);
    int calculateInsertionCost1(const std::vector<Route>& route, int customerId);
    Solution buildRoutes(int period, int dmax, std::vector<std::vector<int>>& currentInventory);
    Solution localSearch(Solution& solution);
    Solution findBestSolution(int n, int dmax);
    double calculateRouteCost(const std::vector<Route>& route);
    void printSeed();
};

#endif // SOLVER_H
