#include "Solver.h"
#include <algorithm>
#include <iostream>

Solver::Solver(const IRP& irp) : irp(irp) {
    std::random_device rd;
    rng.seed(rd());
    printSeed();
}

void Solver::printSeed() {
    std::cout << "Random seed: " << rng() << std::endl;
}

void Solver::updateInventory(int period, int customerId, std::vector<std::vector<int>>& currentInventory, int deliveryAmount) {
    std::cout << "Updating inventory for Customer " << customerId << " at Period " << period << std::endl;
    if (period < irp.customers[customerId].demand.size() && period < currentInventory[customerId].size()) {
        int newInventory = currentInventory[customerId][period] - irp.customers[customerId].demand[period] + deliveryAmount;
        currentInventory[customerId][period] = newInventory;
    } else {
        std::cerr << "Error: Period " << period << " is out of range for customer demand or inventory." << std::endl;
    }
}

int Solver::calculateInsertionCost1(const std::vector<Route>& route, int customerId) {
    double minDistance = std::numeric_limits<double>::max();
    int bestPosition = -1;
    for (int i = 0; i < route.size(); ++i) {
        int lastCustomer = route[i].route.first;
        double distance = irp.costMatrix[lastCustomer][customerId];
        if (distance < minDistance) {
            minDistance = distance;
            bestPosition = i;
        }
    }
    return bestPosition;
}

double Solver::calculateRouteCost(const std::vector<Route>& route) {
    double cost = 0.0;
    for (size_t i = 0; i < route.size() - 1; ++i) {
        cost += irp.costMatrix[route[i].route.first][route[i + 1].route.first];
    }
    return cost;
}

Solution Solver::buildRoutes(int period, int dmax, std::vector<std::vector<int>>& currentInventory) {
    Solution solution(irp);
    solution.currentInventory = currentInventory;

    // Reset vehicle capacities and routes at the beginning of each period
    std::vector<Vehicle> vehicles;
    for (int i = 0; i < irp.Fleet_Size; ++i) {
        vehicles.emplace_back(i, irp.Capacity, irp.depots[0].id);
        solution.vehicleRoutes[period].emplace_back(); // Add a new vector<Route> for each vehicle
        // Start and end each vehicle's route at the depot
        solution.vehicleRoutes[period][i].push_back(Route(irp.depots[0].id, 0, 0.0));
    }

    std::vector<int> candidateList;
    std::vector<bool> customerServed(irp.customers.size(), false); // Track if a customer is already served

    for (int i = 0; i < irp.customers.size(); ++i) {
        if (currentInventory[i][period] - irp.customers[i].demand[period] < 0) { // Needs delivery if inventory - demand < 0
            candidateList.push_back(irp.customers[i].id);
        }
    }

    // Shuffle candidate list to add randomness
    std::shuffle(candidateList.begin(), candidateList.end(), rng);

    while (!candidateList.empty()) {
        // Find the nearest vehicle with sufficient capacity
        int nearestVehicleIndex = -1;
        int nearestCustomerIndex = -1;
        double minDistance = std::numeric_limits<double>::max();

        for (int j = 0; j < candidateList.size(); ++j) {
            int customerId = candidateList[j];
            if (customerId >= irp.customers.size() || period >= irp.customers[customerId].demand.size() || customerServed[customerId]) {
                std::cerr << "Error: Invalid customer or period index or customer already served." << std::endl;
                continue;
            }

            for (int i = 0; i < vehicles.size(); ++i) {
                double distance = irp.costMatrix[vehicles[i].currentLocation][customerId];
                if (distance < minDistance && vehicles[i].currentCapacity >= irp.customers[customerId].demand[period]) {
                    minDistance = distance;
                    nearestVehicleIndex = i;
                    nearestCustomerIndex = j;
                }
            }

            if (nearestVehicleIndex != -1) break;
        }

        if (nearestVehicleIndex == -1 || nearestCustomerIndex == -1) {
            std::cerr << "No feasible routes found. Exiting loop." << std::endl;
            break; // No more feasible routes
        }

        // Calculate random delivery amount
        int customerId = candidateList[nearestCustomerIndex];
        int d0 = irp.customers[customerId].demand[period];
        int dmaxAdjusted = 0;
        for (int k = 0; k < dmax && (period + k) < irp.nPeriods; ++k) {
            dmaxAdjusted += irp.customers[customerId].demand[period + k];
        }
        std::uniform_int_distribution<int> dist(0, dmaxAdjusted - d0);
        int randDelivery = dist(rng);
        int deliveryAmount = std::min(vehicles[nearestVehicleIndex].currentCapacity, d0 + randDelivery);
        deliveryAmount = std::min(deliveryAmount, irp.customers[customerId].maxLevelInv - currentInventory[customerId][period] + irp.customers[customerId].demand[period]);

        // Visit the nearest customer
        vehicles[nearestVehicleIndex].visitCustomer(customerId, deliveryAmount);
        double routeCost = irp.costMatrix[vehicles[nearestVehicleIndex].currentLocation][customerId];
        updateInventory(period, customerId, currentInventory, deliveryAmount);
        candidateList.erase(candidateList.begin() + nearestCustomerIndex);
        customerServed[customerId] = true; // Mark customer as served

        // Update the current period cost
        vehicles[nearestVehicleIndex].currentLocation = customerId;

        // Record the route
        solution.vehicleRoutes[period][nearestVehicleIndex].push_back(Route(customerId, deliveryAmount, routeCost));
    }

    // Ensure all vehicles return to the depot
    for (auto& vehicle : vehicles) {
        double routeCost = irp.costMatrix[vehicle.currentLocation][irp.depots[0].id];
        vehicle.returnToDepot(irp.depots[0].id);
        solution.vehicleRoutes[period][vehicle.id].push_back(Route(irp.depots[0].id, 0, routeCost));
    }

    solution.currentInventory = currentInventory;
    solution.calculateCosts();
    return solution;
}

Solution Solver::solve(int dmax) {
    std::vector<std::vector<int>> currentInventory(irp.nCustomers, std::vector<int>(irp.nPeriods, 0));
    for (int i = 0; i < irp.nCustomers; ++i) {
        currentInventory[i][0] = irp.customers[i].initialInv;
    }

    Solution bestSolution(irp);
    for (int t = 0; t < irp.nPeriods; ++t) {
        std::cout << "Building routes for Period " << t << std::endl;
        Solution periodSolution = buildRoutes(t, dmax, currentInventory);
        bestSolution.vehicleRoutes[t] = periodSolution.vehicleRoutes[t];
        bestSolution.currentInventory = periodSolution.currentInventory;
    }

    bestSolution.calculateCosts();
    return bestSolution;
}

Solution Solver::localSearch(Solution& solution) {
    for (int t = 0; t < irp.nPeriods; ++t) {
        for (auto& vehicle : solution.vehicleRoutes[t]) {
            int n = vehicle.size();
            if (n <= 2) continue; // No need to optimize routes with less than 3 points

            bool improvement = true;
            while (improvement) {
                improvement = false;
                for (int i = 1; i < n - 2; ++i) {  // Skip the first and last element (depot)
                    for (int j = i + 1; j < n - 1; ++j) {  // Skip the first and last element (depot)
                        if (j - i == 1) continue; // Skip adjacent points

                        // Calculate the cost of the current and new segments
                        double oldCost = irp.costMatrix[vehicle[i - 1].route.first][vehicle[i].route.first] +
                                         irp.costMatrix[vehicle[j].route.first][vehicle[j + 1].route.first];

                        double newCost = irp.costMatrix[vehicle[i - 1].route.first][vehicle[j].route.first] +
                                         irp.costMatrix[vehicle[i].route.first][vehicle[j + 1].route.first];

                        if (newCost < oldCost) {
                            std::reverse(vehicle.begin() + i, vehicle.begin() + j + 1);
                            improvement = true;

                            // Recalculate the route cost
                            double routeCost = calculateRouteCost(vehicle);
                            for (auto& route : vehicle) {
                                route.routeCost = routeCost;
                            }
                        }
                    }
                }
            }
        }
    }

    solution.calculateCosts();
    return solution;
}

Solution Solver::findBestSolution(int n, int dmax) {
    Solution bestSolution(irp);
    double bestCost = std::numeric_limits<double>::max();

    for (int i = 0; i < n; ++i) {
        Solution solution = solve(dmax);
        solution = localSearch(solution);

        double currentCost = solution.routeCost + solution.inventoryCost;
        if (currentCost < bestCost) {
            bestCost = currentCost;
            bestSolution = solution;
        }
    }

    return bestSolution;
}
