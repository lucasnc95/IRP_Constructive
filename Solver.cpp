#include "Solver.h"
#include <algorithm>
#include <random>
#include <vector>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>

Solver::Solver(const IRP& irp) : irp(irp) {
    std::random_device rd;
    rng.seed(rd());
}

void Solver::updateInventory(int period, int customerId, std::vector<std::vector<int>>& currentInventory, int deliveryAmount) {
    if (customerId < 1 || customerId > irp.nCustomers || period >= irp.nPeriods) {
        return;
    }
    int newInventory = currentInventory[customerId][period] - irp.customers[customerId-1].demand[period] + deliveryAmount;
    currentInventory[customerId][period] = newInventory;
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
    double alpha = 0.4;
    // Reset vehicle capacities and routes at the beginning of each period
    std::vector<Vehicle> vehicles;
    for (int i = 0; i < irp.Fleet_Size; ++i) {
        vehicles.emplace_back(i, irp.Capacity, irp.depots[0].id);
        solution.vehicleRoutes[period].emplace_back(); // Add a new Route for each vehicle
        // Start each vehicle's route at the depot
        solution.vehicleRoutes[period][i].route.push_back({irp.depots[0].id, 0});
    }

    std::vector<int> candidateList;
    std::vector<bool> customerServed(irp.nCustomers + 1, false); // Track if a customer is already served

    // Build candidate list based on demand and alpha probability
    for (int i = 1; i <= irp.nCustomers; ++i) {
        if (currentInventory[i][period] - irp.customers[i-1].demand[period] < 0) { // Needs delivery if inventory - demand < 0
            candidateList.push_back(irp.customers[i-1].id);
        } else {
            // Check probability to add the customer even if inventory is sufficient
            if (alpha > 0.0 && alpha < 1.0) {
                double randomValue = static_cast<double>(rand()) / RAND_MAX;
                if (randomValue < alpha) {
                    candidateList.push_back(irp.customers[i-1].id);
                }
            }
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
            if (customerId < 1 || customerId > irp.nCustomers || period >= irp.nPeriods || customerServed[customerId]) {
                continue;
            }

            for (int i = 0; i < vehicles.size(); ++i) {
                double distance = irp.costMatrix[vehicles[i].currentLocation][customerId];
                if (distance < minDistance && vehicles[i].currentCapacity >= irp.customers[customerId-1].demand[period]) {
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
        int d0 = irp.customers[customerId-1].demand[period];
        int dmaxAdjusted = 0;
        for (int k = 0; k < dmax && (period + k) < irp.nPeriods; ++k) {
            dmaxAdjusted += irp.customers[customerId-1].demand[period + k];
        }
        std::uniform_int_distribution<int> dist(0, dmaxAdjusted - d0);
        int randDelivery = dist(rng);
        int deliveryAmount = std::min(vehicles[nearestVehicleIndex].currentCapacity, d0 + randDelivery);
        deliveryAmount = std::min(deliveryAmount, irp.customers[customerId-1].maxLevelInv - currentInventory[customerId][period] + irp.customers[customerId-1].demand[period]);

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
    std::vector<std::vector<int>> currentInventory(irp.nCustomers + 1, std::vector<int>(irp.nPeriods, 0));
    for (int i = 1; i <= irp.nCustomers; ++i) {
        currentInventory[i][0] = irp.customers[i-1].initialInv;
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
        double bestCost = route.routeCost;

        for (int i = 1; i < n - 2; ++i) {  // Skip the first and last element (depot)
            for (int j = i + 1; j < n - 1; ++j) {  // Skip the first and last element (depot)
                if (j - i == 1) continue; // Skip adjacent points

                // Calculate the cost of the current and new segments
                double oldCost = irp.costMatrix[route.route[i - 1].first][route.route[i].first] +
                                 irp.costMatrix[route.route[j].first][route.route[j + 1].first];

                double newCost = irp.costMatrix[route.route[i - 1].first][route.route[j].first] +
                                 irp.costMatrix[route.route[i].first][route.route[j + 1].first];

                double costImprovement = oldCost - newCost;

                if (costImprovement > 0) {
                    std::reverse(route.route.begin() + i, route.route.begin() + j + 1);
                    improvement = true;

                    // Recalculate the route cost
                    route.routeCost = calculateRouteCost(route);
                    bestCost = route.routeCost;
                }
            }
        }

        if (!improvement) {
            route.routeCost = bestCost;
        }
    }
}

void Solver::swap(Route& route) {
    int n = route.route.size();
    if (n <= 2) return; // No need to optimize routes with less than 3 points

    bool improvement = true;
    while (improvement) {
        improvement = false;
        double bestCost = route.routeCost;

        for (int i = 1; i < n - 1; ++i) {  // Skip the first and last element (depot)
            for (int j = i + 1; j < n - 1; ++j) {  // Skip the first and last element (depot)
                // Swap customers
                std::swap(route.route[i], route.route[j]);

                // Recalculate the route cost
                double newCost = calculateRouteCost(route);
                double costImprovement = route.routeCost - newCost;

                if (costImprovement > 0) {
                    route.routeCost = newCost;
                    improvement = true;
                } else {
                    // Swap back if no improvement
                    std::swap(route.route[i], route.route[j]);
                }
            }
        }

        if (!improvement) {
            route.routeCost = bestCost;
        }
    }
}

void Solver::relocate(Route& route) {
    int n = route.route.size();
    if (n <= 2) return; // No need to optimize routes with less than 3 points

    bool improvement = true;
    while (improvement) {
        improvement = false;
        double bestCost = route.routeCost;

        for (int i = 1; i < n - 1; ++i) {  // Skip the first and last element (depot)
            auto customer = route.route[i];
            for (int j = 1; j < n - 1; ++j) {
                if (i == j) continue;

                // Remove and insert customer at new position
                route.route.erase(route.route.begin() + i);
                route.route.insert(route.route.begin() + j, customer);

                // Recalculate the route cost
                double newCost = calculateRouteCost(route);
                double costImprovement = route.routeCost - newCost;

                if (costImprovement > 0) {
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

        if (!improvement) {
            route.routeCost = bestCost;
        }
    }
}

void Solver::exchangeRoutes(Route& route1, Route& route2) {
    bool improvement = false;
    double bestImprovement = 0.0;
    int best_i = -1, best_j = -1;

    for (int i = 1; i < route1.route.size() - 1; ++i) {
        for (int j = 1; j < route2.route.size() - 1; ++j) {
            auto& customer1 = route1.route[i];
            auto& customer2 = route2.route[j];

            // Check if exchange is feasible
            if (route1.remainingCapacity + customer1.second - customer2.second >= 0 &&
                route2.remainingCapacity + customer2.second - customer1.second >= 0) {

                // Calculate potential new costs
                std::swap(customer1, customer2);
                double newCost1 = calculateRouteCost(route1);
                double newCost2 = calculateRouteCost(route2);
                std::swap(customer1, customer2);

                double costImprovement = (route1.routeCost + route2.routeCost) - (newCost1 + newCost2);

                if (costImprovement > bestImprovement) {
                    bestImprovement = costImprovement;
                    best_i = i;
                    best_j = j;
                    improvement = true;
                }
            }
        }
    }

    if (improvement) {
        // Apply the best exchange found
        std::swap(route1.route[best_i], route2.route[best_j]);
        route1.routeCost = calculateRouteCost(route1);
        route2.routeCost = calculateRouteCost(route2);
        route1.remainingCapacity += route1.route[best_i].second - route2.route[best_j].second;
        route2.remainingCapacity += route2.route[best_j].second - route1.route[best_i].second;
    }
}

void Solver::relocateRoutes(Route& route1, Route& route2) {
    bool improvement = false;
    double bestImprovement = 0.0;
    int best_i = -1, best_j = -1;

    for (int i = 1; i < route1.route.size() - 1; ++i) {
        auto customer = route1.route[i];

        for (int j = 1; j < route2.route.size(); ++j) {
            if (route2.remainingCapacity >= customer.second) {
                // Calculate potential new costs
                route1.route.erase(route1.route.begin() + i);
                route2.route.insert(route2.route.begin() + j, customer);
                double newCost1 = calculateRouteCost(route1);
                double newCost2 = calculateRouteCost(route2);
                route2.route.erase(route2.route.begin() + j);
                route1.route.insert(route1.route.begin() + i, customer);

                double costImprovement = (route1.routeCost + route2.routeCost) - (newCost1 + newCost2);

                if (costImprovement > bestImprovement) {
                    bestImprovement = costImprovement;
                    best_i = i;
                    best_j = j;
                    improvement = true;
                }
            }
        }
    }

    if (improvement) {
        // Apply the best relocation found
        auto customer = route1.route[best_i];
        route1.route.erase(route1.route.begin() + best_i);
        route2.route.insert(route2.route.begin() + best_j, customer);
        route1.routeCost = calculateRouteCost(route1);
        route2.routeCost = calculateRouteCost(route2);
        route1.remainingCapacity += customer.second;
        route2.remainingCapacity -= customer.second;
    }
}

void Solver::swapRoutes(Route& route1, Route& route2) {
    bool improvement = false;
    double bestImprovement = 0.0;
    int best_i = -1, best_j = -1;

    for (int i = 1; i < route1.route.size() - 1; ++i) {
        for (int j = 1; j < route2.route.size() - 1; ++j) {
            auto& customer1 = route1.route[i];
            auto& customer2 = route2.route[j];

            // Check if swap is feasible
            if (route1.remainingCapacity + customer1.second - customer2.second >= 0 &&
                route2.remainingCapacity + customer2.second - customer1.second >= 0) {

                // Calculate potential new costs
                std::swap(customer1, customer2);
                double newCost1 = calculateRouteCost(route1);
                double newCost2 = calculateRouteCost(route2);
                std::swap(customer1, customer2);

                double costImprovement = (route1.routeCost + route2.routeCost) - (newCost1 + newCost2);

                if (costImprovement > bestImprovement) {
                    bestImprovement = costImprovement;
                    best_i = i;
                    best_j = j;
                    improvement = true;
                }
            }
        }
    }

    if (improvement) {
        // Apply the best swap found
        std::swap(route1.route[best_i], route2.route[best_j]);
        route1.routeCost = calculateRouteCost(route1);
        route2.routeCost = calculateRouteCost(route2);
        route1.remainingCapacity += route1.route[best_i].second - route2.route[best_j].second;
        route2.remainingCapacity += route2.route[best_j].second - route1.route[best_i].second;
    }
}

Solution Solver::localSearchAcrossRoutes(Solution& solution, int iterations) {
    std::random_device rd;
    std::mt19937 g(rd());

    for (int n = 0; n < iterations; ++n) {
        for (int t = 0; t < irp.nPeriods; ++t) {
            std::shuffle(solution.vehicleRoutes[t].begin(), solution.vehicleRoutes[t].end(), g);

            for (size_t i = 0; i < solution.vehicleRoutes[t].size(); ++i) {
                for (size_t j = i + 1; j < solution.vehicleRoutes[t].size(); ++j) {
                    int localSearchType = std::uniform_int_distribution<int>(1, 3)(g);

                    if (localSearchType == 1) {
                        exchangeRoutes(solution.vehicleRoutes[t][i], solution.vehicleRoutes[t][j]);
                    } else if (localSearchType == 2) {
                        relocateRoutes(solution.vehicleRoutes[t][i], solution.vehicleRoutes[t][j]);
                    } else if (localSearchType == 3) {
                        swapRoutes(solution.vehicleRoutes[t][i], solution.vehicleRoutes[t][j]);
                    }
                }
            }
        }
    }

    solution.calculateCosts();
    return solution;
}

Solution Solver::localSearch(Solution& solution, int iterations) {
    std::random_device rd;
    std::mt19937 g(rd());

    for (int n = 0; n < iterations; ++n) {
        for (int t = 0; t < irp.nPeriods; ++t) {
            for (auto& vehicle : solution.vehicleRoutes[t]) {
                double bestCost = vehicle.routeCost;

                // Randomly apply twoOpt, swap, and relocate with best improvement criteria
                int localSearchType = std::uniform_int_distribution<int>(1, 3)(g);

                if (localSearchType == 1) {
                    twoOpt(vehicle);
                } else if (localSearchType == 2) {
                    swap(vehicle);
                } else if (localSearchType == 3) {
                    relocate(vehicle);
                }

                // Check if there is any improvement
                if (vehicle.routeCost < bestCost) {
                    bestCost = vehicle.routeCost;
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
        solution = localSearch(solution, 100);
        solution = localSearchAcrossRoutes(solution, 100);

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

    // Print all customers and depots with their coordinates
    file << "Depots:\n";
    for (const auto& depot : irp.depots) {
        file << depot.id << " " << depot.x << " " << depot.y << "\n";
    }

    file << "Customers:\n";
    for (const auto& customer : irp.customers) {
        file << customer.id << " " << customer.x << " " << customer.y << "\n";
    }

    // Print the routes for each period and each vehicle
    for (int t = 0; t < irp.nPeriods; ++t) {
        file << "Period " << t << ":\n";
        for (const auto& vehicle : solution.vehicleRoutes[t]) {
            file << "  Vehicle " << (&vehicle - &solution.vehicleRoutes[t][0]) << ": ";
            for (const auto& route : vehicle.route) {
                file << route.first << " ";
            }
            file << "\n";
        }
    }

    // Print the inventory levels for each period and each customer
    file << "Inventory levels:\n";
    for (int t = 0; t < irp.nPeriods; ++t) {
        file << "Period " << t << ":\n";
        for (int i = 1; i <= irp.nCustomers; ++i) {
            file << "  Customer " << irp.customers[i-1].id << ": " << solution.currentInventory[i][t] << "\n";
        }
    }

    file.close();

    // Call the Python script to plot the solution
    std::string command = "python3 plot_solution.py " + filename;
    system(command.c_str());
}
