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
    Solution solve(int dmax, float alpha);
    Solution buildRoutes(int period, int dmax, std::vector<std::vector<int>>& currentInventory, float alpha);
    Solution buildRoutesEnhanced(int period, int dmax, std::vector<std::vector<int>>& currentInventory, float alpha);
    void updateInventory(int period, int customerId, std::vector<std::vector<int>>& currentInventory, int deliveryAmount);
    double calculateRouteCost(const Route& route);
    Solution localSearch(Solution& solution, int iterations);
    Solution localSearchAcrossRoutes(Solution& solution, int iterations);
    Solution localSearchRandomized(Solution& solution, int maxNoImprovement,unsigned seed);
    void executeLocalSearch(int searchIndex, Solution& solution);
    Solution findBestSolution(int n, int maxDmax, unsigned seed, int  maxNoImp);
    void saveSolutionToFile(const Solution& solution, const std::string& filename);
    double calculateInsertionCost(const Route& route, size_t pos, int customerId) const;
    void performInsertion(Route& route, size_t pos, int customerId, int quantity);
    bool validateRouteConsistency(const Route& route) const;
    double calculateCandidateScore(int customerId, int period, const std::vector<std::vector<int>>& inventory) const;
    void reconstructRoute(Route& route) const;
    size_t findBestInsertionPosition(const Route& route, int customerId) const;   
    double calculateRouteCost(const Route& route) const {
        double totalCost = 0.0;
        
        // Calcula o custo entre cada par consecutivo de nós na rota
        for (size_t i = 1; i < route.route.size(); ++i) {
            int from = route.route[i-1].first;
            int to = route.route[i].first;
            totalCost += irp.costMatrix[from][to];
        }
        
        return totalCost;
    }

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
    bool supplyInsertion(Solution& solution, int periodFrom, int periodTo, int customerId, int quantity);
    bool supplyRemoval(Solution& solution, int period, int customerId, int quantity);
    bool shiftDelivery(Solution& solution, int periodFrom, int periodTo, int customerId);
    bool localSearchInventory(Solution& solution);
    bool validateSolution(Solution& solution);
    bool validateInventory(Solution& solution);
    void initializeNewVehicleRoute(Solution& solution, int period, int customerId, int quantity);


private:
    const IRP& irp;
    std::mt19937 rng;
};


#endif // SOLVER_H
