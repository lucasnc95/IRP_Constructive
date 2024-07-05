#ifndef SOLUTION_H
#define SOLUTION_H

#include "IRP.h"
#include <vector>
#include <utility>
//#include <Route.h>

class Route {
public:
    std::pair<int, int> route; // (customer_id, delivery_amount)
    double routeCost;

    Route(int customerId, int deliveryAmount, double cost)
        : route(customerId, deliveryAmount), routeCost(cost) {}
};

class Solution {
public:
    const IRP& irp;
    std::vector<std::vector<std::vector<Route>>> vehicleRoutes;
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
