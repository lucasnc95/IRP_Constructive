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
            double vehicleRouteCost = 0.0;
            for (size_t i = 0; i < vehicle.route.size() - 1; ++i) {
                vehicleRouteCost += irp.costMatrix[vehicle.route[i].first][vehicle.route[i + 1].first];
            }
            routeCost += vehicleRouteCost;
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
    for (int i = 1; i <= irp.nCustomers; ++i) {  // Clientes vão de 1 a irp.nCustomers
        for (int t = 0; t < irp.nPeriods; ++t) {
            if (currentInventory[i][t] < 0) {
                std::cerr << "Customer " << i << " has negative inventory at period " << t << std::endl;
                return false;
            }
            if (currentInventory[i][t] > irp.customers[i-1].maxLevelInv) {
                std::cerr << "Customer " << i << " exceeds max inventory level at period " << t << std::endl;
                return false;
            }
        }
    }

    // Check for multiple visits to the same customer in the same period
    for (int t = 0; t < irp.nPeriods; ++t) {
        std::vector<int> visits(irp.nCustomers + 1, 0);  // Clientes vão de 1 a irp.nCustomers
        for (const auto& vehicle : vehicleRoutes[t]) {
            int vehicleLoad = 0;
            if (vehicle.route.front().first != 0 || vehicle.route.back().first != 0) {
                std::cerr << "Vehicle does not start and end at the depot in period " << t << std::endl;
                return false;
            }
            for (const auto& route : vehicle.route) {
                int customerId = route.first;
                if (customerId != 0) {  // Only depot 0 can be visited more than once
                    visits[customerId]++;
                    if (visits[customerId] > 1) {
                        std::cerr << "Customer " << customerId << " visited more than once in period " << t << std::endl;
                        return false;
                    }
                }
                vehicleLoad += route.second;
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
            for (const auto& route : vehicle.route) {
                std::cout << "(ID " << route.first << ",  " << route.second << ") ";
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


// void Route::addDelivery(int customerId, int quantity) {
//     bool customerFound = false;

//     for (auto& stop : route) {
//         if (stop.first == customerId) {
//             stop.second += quantity;
//             customerFound = true;
//             break;
//         }
//     }

//     if (!customerFound) {
//         route.push_back({customerId, quantity});
//     }

//     cargaTotal += quantity;
//     remainingCapacity -= quantity;
// }

// // Remove uma quantidade específica de entrega de um cliente
// void Route::removeDelivery(int customerId, int quantity) {
//     for (auto it = route.begin(); it != route.end(); ++it) {
//         if (it->first == customerId) {
//             if (it->second > quantity) {
//                 // Reduz apenas parte da entrega
//                 it->second -= quantity;
//                 cargaTotal -= quantity;
//                 remainingCapacity += quantity;
//             } else {
//                 // Remove completamente a entrega se a quantidade for igual ou maior
//                 cargaTotal -= it->second;
//                 remainingCapacity += it->second;
//                 route.erase(it);
//             }
//             return; // A operação foi concluída
//         }
//     }

//   //  std::cerr << "Erro: Cliente " << customerId << " não encontrado na rota." << std::endl;
// }


double Solution::getTotalCost(){


calculateCosts();


return routeCost + inventoryCost;

}



void Route::addDelivery(int customerId, int quantity) {
    route.insert(route.end() - 1, {customerId, quantity});
    cargaTotal += quantity;
    remainingCapacity -= quantity;
}


void Route::removeDelivery(size_t position) {
    if (position < route.size()) {
        remainingCapacity += route[position].second;
        cargaTotal -= route[position].second;
        route.erase(route.begin() + position);
    }
}