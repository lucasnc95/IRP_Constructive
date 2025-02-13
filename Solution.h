#ifndef SOLUTION_H
#define SOLUTION_H

#include "IRP.h"
#include <vector>
#include <utility>

class Route {
public:
    std::vector<std::pair<int, int>> route; // <customer ID, delivered quantity>
    int cargaTotal;
    double routeCost;
    int remainingCapacity; // New variable to store remaining capacity

    void printRoute() const {
        std::cout << "Route details:" << std::endl;
        for (const auto& stop : route) {
            std::cout << "(ID " << stop.first << ",  " << stop.second << ") ";
        }
       std::cout<<std::endl;
    }

    Route() : cargaTotal(0), routeCost(0.0), remainingCapacity(0) {}
    void addDelivery(int customerId, int quantity);
    void removeDelivery(int customerId, int quantity);
};


class Solution {
public:
    const IRP& irp;
    std::vector<std::vector<Route>> vehicleRoutes;
    std::vector<std::vector<int>> currentInventory;
    double routeCost;
    double inventoryCost;
    double bestAlpha = -1.0;
    int bestDmax = -1;
    Solution(const IRP& irp);
    Solution(const Solution& other);
    Solution& operator=(const Solution& other);

    void calculateCosts();
    bool isFeasible() const;
    void printSolution() const;
    void printInventoryLevels(int period) const;
    double getTotalCost();
    
};

#endif // SOLUTION_H

