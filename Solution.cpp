#include "Solution.h"
#include <iostream>

Solution::Solution(const IRP& irp) : irp(irp), routeCost(0), inventoryCost(0) {
    vehicleRoutes.resize(irp.nPeriods);
    currentInventory.resize(irp.nCustomers, std::vector<int>(irp.nPeriods, 0));
    for (int i = 0; i < irp.nCustomers; ++i) {
        currentInventory[i][0] = irp.customers[i].initialInv;
    }
}

Solution::Solution(const Solution& other) : irp(other.irp), vehicleRoutes(other.vehicleRoutes), currentInventory(other.currentInventory), routeCost(other.routeCost), inventoryCost(other.inventoryCost) {
    // Copy constructor implementation
}

Solution& Solution::operator=(const Solution& other) {
    if (this == &other) {
        return *this; // Handle self-assignment
    }
    vehicleRoutes = other.vehicleRoutes;
    currentInventory = other.currentInventory;
    routeCost = other.routeCost;
    inventoryCost = other.inventoryCost;
    return *this;
}

void Solution::calculateCosts() {
    routeCost = 0.0;
    inventoryCost = 0.0;

    // Calculate travel costs
    for (int t = 0; t < irp.nPeriods; ++t) {
        for (const auto& vehicle : vehicleRoutes[t]) {
            for (const auto& route : vehicle) {
                routeCost += route.routeCost;
            }
        }
    }

    // Calculate inventory holding costs
    for (int t = 0; t < irp.nPeriods; ++t) {
        for (int i = 0; i < irp.customers.size(); ++i) {
            if (t < currentInventory[i].size()) {
                inventoryCost += currentInventory[i][t] * irp.customers[i].invCost;
            }
        }
    }
}

bool Solution::isFeasible() const {
    // Check inventory levels and customer visits
    for (int i = 0; i < irp.customers.size(); ++i) {
        for (int t = 0; t < irp.nPeriods; ++t) {
            if (currentInventory[i][t] < 0) {
                std::cerr << "Customer " << i << " has negative inventory at period " << t << std::endl;
                return false;
            }
            if (currentInventory[i][t] > irp.customers[i].maxLevelInv) {
                std::cerr << "Customer " << i << " exceeds max inventory level at period " << t << std::endl;
                return false;
            }
        }
    }

    // Check for multiple visits to the same customer in the same period
    for (int t = 0; t < irp.nPeriods; ++t) {
        std::vector<int> visits(irp.customers.size(), 0);
        for (const auto& vehicle : vehicleRoutes[t]) {
            int vehicleLoad = 0;
            if (vehicle.front().route.first != 0 || vehicle.back().route.first != 0) {
                std::cerr << "Vehicle does not start and end at the depot in period " << t << std::endl;
                return false;
            }
            for (const auto& route : vehicle) {
                int customerId = route.route.first;
                if (customerId != 0) {  // Only depot 0 can be visited more than once
                    visits[customerId]++;
                    if (visits[customerId] > 1) {
                        std::cerr << "Customer " << customerId << " visited more than once in period " << t << std::endl;
                        return false;
                    }
                }
                vehicleLoad += route.route.second;
            }
            // Check if vehicle load exceeds capacity
            if (vehicleLoad > irp.Capacity) {
                std::cerr << "Vehicle load exceeds capacity in period " << t << std::endl;
                return false;
            }
        }
    }

    return true;
}

void Solution::printSolution() const {
    for (int t = 0; t < irp.nPeriods; ++t) {
        std::cout << "Period " << t << ":\n";
        for (const auto& vehicle : vehicleRoutes[t]) {
            std::cout << "  Route: ";
            for (const auto& route : vehicle) {
                std::cout << "<Cliente " << route.route.first << ", Quantidade entregue: " << route.route.second << "> ";
            }
            std::cout << "\n";
        }
    }
}

void Solution::printInventoryLevels(int period) const {
    std::cout << "Inventory levels at period " << period << ":\n";
    for (int i = 0; i < irp.customers.size(); ++i) {
        std::cout << "  Customer " << irp.customers[i].id << ": " << currentInventory[i][period] << "\n";
    }
}
