#ifndef SOLUTION_H
#define SOLUTION_H

#include "IRP.h"
#include <vector>
#include <utility>

class Route {
public:
    std::vector<std::pair<int, int>> route; // Vector of (customer_id, delivery_amount) pairs
    double routeCost;

    Route() : routeCost(0.0) {}
};

class Solution {
public:
    const IRP& irp;
    std::vector<std::vector<Route>> vehicleRoutes;
    std::vector<std::vector<int>> currentInventory;
    double routeCost;
    double inventoryCost;

    Solution(const IRP& irp);
    Solution(const Solution& other);
    Solution& operator=(const Solution& other);

    void calculateCosts();
    bool isFeasible() const;
    void printSolution() const;
    void printInventoryLevels(int period) const;
};

#endif // SOLUTION_H
