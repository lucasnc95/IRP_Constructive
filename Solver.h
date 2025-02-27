#ifndef SOLVER_H
#define SOLVER_H

#include "IRP.h"
#include "Solution.h"
#include "Vehicle.h"
#include <random>
class Solver {
public:
    Solver(const IRP& irp);
    Solution solve(int dmax, float alpha);
    Solution buildRoutes(int period, int dmax, std::vector<std::vector<int>>& currentInventory, float alpha);    
    void adjustInventoryDistribution(Solution& solution, int period, std::vector<std::vector<int>>& currentInventory);
    void allocateInventoryWithConstraints(Solution& solution, std::vector<std::vector<int>>& inventory, int dmax);
    void updateInventory(int period, int customerId, std::vector<std::vector<int>>& currentInventory, int deliveryAmount);
    double calculateRouteCost(const Route& route);
    Solution localSearch(Solution& solution, int iterations);
    Solution localSearchAcrossRoutes(Solution& solution, int iterations);
    Solution localSearchRandomized(Solution& solution, int maxNoImprovement,unsigned seed);
    void executeLocalSearch(int searchIndex, Solution& solution);
    Solution findBestSolution(int n, int maxDmax, unsigned seed, int  maxNoImp);
    void saveSolutionToFile(const Solution& solution, const std::string& filename);
    void allocateInventory(Solution& solution, std::vector<std::vector<int>>& inventory,int dmax);
    void twoOpt(Route& route);
    void exchangeOneOne(Route& route);
    void exchangeTwoOne(Route& route);
    void exchangeTwoTwo(Route& route);
    void reinsertion(Route& route);
    void swapOneOne(Route& route1, Route& route2);
    void shiftOneZero(Route& route1, Route& route2);
    void swapRoutes(Route& route1, Route& route2);
    void swapTwoTwoInterRoute(Route& route1, Route& route2);
    void swapTwoOneInterRoute(Route& route1, Route& route2);
    void shiftTwoZeroInterRoute(Route& route1, Route& route2);
    bool localSearchInventory(Solution& solution);
 


private:
    const IRP& irp;
    std::mt19937 rng;
};


#endif // SOLVER_H
