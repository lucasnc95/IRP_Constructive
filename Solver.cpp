#include "Solver.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>

Solver::Solver(const IRP& irp) : irp(irp) {
    std::random_device rd;
    rng.seed(rd());
}

void Solver::updateInventory(int period, int customerId, std::vector<std::vector<int>>& currentInventory, int deliveryAmount) {
    if (period < irp.customers[customerId].demand.size() && period < currentInventory[customerId].size()) {
        int newInventory = currentInventory[customerId][period] - irp.customers[customerId].demand[period] + deliveryAmount;
        currentInventory[customerId][period] = newInventory;
    }
}

int Solver::calculateInsertionCost1(const std::vector<Route>& route, int customerId) {
    double minDistance = std::numeric_limits<double>::max();
    int bestPosition = -1;
    for (int i = 0; i < route.size(); ++i) {
        int lastCustomer = route[i].route.back().first;
        double distance = irp.costMatrix[lastCustomer][customerId];
        if (distance < minDistance) {
            minDistance = distance;
            bestPosition = i;
        }
    }
    return bestPosition;
}

double Solver::calculateRouteCost(const Route& route) {
    double cost = 0.0;
    for (size_t i = 0; i < route.route.size() - 1; ++i) {
        cost += irp.costMatrix[route.route[i].first][route.route[i + 1].first];
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
        solution.vehicleRoutes[period].emplace_back(); // Add a new Route for each vehicle
        // Start each vehicle's route at the depot
        solution.vehicleRoutes[period][i].route.push_back({irp.depots[0].id, 0});
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
        solution.vehicleRoutes[period][nearestVehicleIndex].route.push_back({customerId, deliveryAmount});
        solution.vehicleRoutes[period][nearestVehicleIndex].routeCost += routeCost;
    }

    // Ensure all vehicles return to the depot
    for (auto& vehicle : vehicles) {
        double routeCost = irp.costMatrix[vehicle.currentLocation][irp.depots[0].id];
        vehicle.returnToDepot(irp.depots[0].id);
        solution.vehicleRoutes[period][vehicle.id].route.push_back({irp.depots[0].id, 0});
        solution.vehicleRoutes[period][vehicle.id].routeCost += routeCost;
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
        Solution periodSolution = buildRoutes(t, dmax, currentInventory);
        bestSolution.vehicleRoutes[t] = periodSolution.vehicleRoutes[t];
        bestSolution.currentInventory = periodSolution.currentInventory;
    }

    bestSolution.calculateCosts();
    return bestSolution;
}

void Solver::twoOpt(Route& route) {
    int n = route.route.size();
    if (n <= 2) return; // No need to optimize routes with less than 3 points

    bool improvement = true;
    while (improvement) {
        improvement = false;
        for (int i = 1; i < n - 2; ++i) {  // Skip the first and last element (depot)
            for (int j = i + 1; j < n - 1; ++j) {  // Skip the first and last element (depot)
                if (j - i == 1) continue; // Skip adjacent points

                // Calculate the cost of the current and new segments
                double oldCost = irp.costMatrix[route.route[i - 1].first][route.route[i].first] +
                                 irp.costMatrix[route.route[j].first][route.route[j + 1].first];

                double newCost = irp.costMatrix[route.route[i - 1].first][route.route[j].first] +
                                 irp.costMatrix[route.route[i].first][route.route[j + 1].first];

                if (newCost < oldCost) {
                    std::reverse(route.route.begin() + i, route.route.begin() + j + 1);
                    improvement = true;

                    // Recalculate the route cost
                    route.routeCost = calculateRouteCost(route);
                }
            }
        }
    }
}

void Solver::swap(Route& route) {
    int n = route.route.size();
    if (n <= 2) return; // No need to optimize routes with less than 3 points

    bool improvement = true;
    while (improvement) {
        improvement = false;
        for (int i = 1; i < n - 1; ++i) {  // Skip the first and last element (depot)
            for (int j = i + 1; j < n - 1; ++j) {  // Skip the first and last element (depot)
                // Swap customers
                std::swap(route.route[i], route.route[j]);

                // Recalculate the route cost
                double newCost = calculateRouteCost(route);
                if (newCost < route.routeCost) {
                    route.routeCost = newCost;
                    improvement = true;
                } else {
                    // Swap back if no improvement
                    std::swap(route.route[i], route.route[j]);
                }
            }
        }
    }
}

void Solver::relocate(Route& route) {
    int n = route.route.size();
    if (n <= 2) return; // No need to optimize routes with less than 3 points

    bool improvement = true;
    while (improvement) {
        improvement = false;
        for (int i = 1; i < n - 1; ++i) {  // Skip the first and last element (depot)
            auto customer = route.route[i];
            for (int j = 1; j < n - 1; ++j) {
                if (i == j) continue;

                // Remove and insert customer at new position
                route.route.erase(route.route.begin() + i);
                route.route.insert(route.route.begin() + j, customer);

                // Recalculate the route cost
                double newCost = calculateRouteCost(route);
                if (newCost < route.routeCost) {
                    route.routeCost = newCost;
                    improvement = true;
                    break;
                } else {
                    // Revert changes if no improvement
                    route.route.erase(route.route.begin() + j);
                    route.route.insert(route.route.begin() + i, customer);
                }
            }
            if (improvement) break;
        }
    }
}

Solution Solver::localSearch(Solution& solution) {
    for (int t = 0; t < irp.nPeriods; ++t) {
        for (auto& vehicle : solution.vehicleRoutes[t]) {
            // Apply 2-opt
            twoOpt(vehicle);

            // Apply swap
            swap(vehicle);

            // Apply relocate
            relocate(vehicle);
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

void Solver::saveSolutionToFile(const Solution& solution, const std::string& filename) {
    std::ofstream file(filename);

    if (!file) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return;
    }

    for (int t = 0; t < irp.nPeriods; ++t) {
        file << "Period " << t << ":\n";
        for (const auto& vehicle : solution.vehicleRoutes[t]) {
            file << "  Vehicle " << &vehicle - &solution.vehicleRoutes[t][0] << ": ";
            for (const auto& route : vehicle.route) {
                file << "(" << route.first << ", " << route.second << ") ";
            }
            file << "\n";
        }
    }

    file << "Inventory levels:\n";
    for (int t = 0; t < irp.nPeriods; ++t) {
        file << "Period " << t << ":\n";
        for (int i = 0; i < irp.customers.size(); ++i) {
            file << "  Customer " << irp.customers[i].id << ": " << solution.currentInventory[i][t] << "\n";
        }
    }

    file.close();

    // Call the Python script to plot the solution
    std::string command = "python3 plot_solution.py " + filename;
    system(command.c_str());
}
