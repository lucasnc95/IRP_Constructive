#ifndef SOLVER_H
#define SOLVER_H

#include "IRP.h"
#include "Solution.h"
#include "Vehicle.h"
//#include <Route.h>
#include <random>
class Solver {
public:
    Solver(const IRP& irp);
    Solution solve(int dmax);
    Solution buildRoutes(int period, int dmax, std::vector<std::vector<int>>& currentInventory);
    void updateInventory(int period, int customerId, std::vector<std::vector<int>>& currentInventory, int deliveryAmount);
    double calculateRouteCost(const Route& route);
    Solution localSearch(Solution& solution, int iterations);
    Solution localSearchAcrossRoutes(Solution& solution, int iterations);
    Solution findBestSolution(int n, int dmax);
    void saveSolutionToFile(const Solution& solution, const std::string& filename);

    void twoOpt(Route& route);
    void swap(Route& route);
    void relocate(Route& route);
    void exchangeRoutes(Route& route1, Route& route2);
    void relocateRoutes(Route& route1, Route& route2);
    void swapRoutes(Route& route1, Route& route2);

private:
    const IRP& irp;
    std::mt19937 rng;
};


#endif // SOLVER_H
