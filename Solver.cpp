#include "Solver.h"
#include <algorithm>
#include <random>
#include <vector>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <functional>
#include <numeric>


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



Solution Solver::solve(int dmax, float alpha) {
    std::vector<std::vector<int>> currentInventory(irp.nCustomers + 1, std::vector<int>(irp.nPeriods, 0));
    for (int i = 1; i <= irp.nCustomers; ++i) {
        currentInventory[i][0] = irp.customers[i-1].initialInv;
    }

    Solution bestSolution(irp);
    for (int t = 0; t < irp.nPeriods; ++t) {
        Solution periodSolution = buildRoutes(t, dmax, currentInventory, alpha);
        bestSolution.vehicleRoutes[t] = periodSolution.vehicleRoutes[t];
        bestSolution.currentInventory = periodSolution.currentInventory;
    }

    bestSolution.calculateCosts();
    return bestSolution;
}


Solution Solver::buildRoutes(int period, int dmax, std::vector<std::vector<int>>& currentInventory, float alpha) {
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
    std::vector<bool> customerServed(irp.nCustomers + 1, false); // Track if a customer is already served

    // Build candidate list based on demand and alpha probability
    for (int i = 1; i <= irp.nCustomers; ++i) {
        if (currentInventory[i][period] - irp.customers[i - 1].demand[period] < 0) { // Needs delivery if inventory - demand < 0
            candidateList.push_back(irp.customers[i - 1].id);
        } else {
            // Add customers randomly based on alpha
            if (alpha > 0.0 && alpha < 1.0) {
                float randomValue = static_cast<float>(rand()) / RAND_MAX;
                if (randomValue < alpha) {
                    candidateList.push_back(irp.customers[i - 1].id);
                }
            }
        }
    }

    // Shuffle candidate list to add randomness
    std::shuffle(candidateList.begin(), candidateList.end(), rng);

    while (!candidateList.empty()) {
        int nearestVehicleIndex = -1;
        int nearestCustomerIndex = -1;
        double minDistance = std::numeric_limits<double>::max();

        // Find nearest vehicle with sufficient capacity
        for (int j = 0; j < candidateList.size(); ++j) {
            int customerId = candidateList[j];
            if (customerId < 1 || customerId > irp.nCustomers || period >= irp.nPeriods || customerServed[customerId]) {
                continue;
            }

            for (int i = 0; i < vehicles.size(); ++i) {
                double distance = irp.costMatrix[vehicles[i].currentLocation][customerId];
                if (distance < minDistance && vehicles[i].currentCapacity >= irp.customers[customerId - 1].demand[period]) {
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

        int customerId = candidateList[nearestCustomerIndex];
        int d0 = irp.customers[customerId - 1].demand[period];
        int dmaxAdjusted = 0;
        for (int k = 0; k < dmax && (period + k) < irp.nPeriods; ++k) {
            dmaxAdjusted += irp.customers[customerId - 1].demand[period + k];
        }

        std::uniform_int_distribution<int> dist(0, dmaxAdjusted - d0);
        int randDelivery = dist(rng);
        int deliveryAmount = std::min(vehicles[nearestVehicleIndex].currentCapacity, d0 + randDelivery);
        deliveryAmount = std::min(deliveryAmount, irp.customers[customerId - 1].maxLevelInv - currentInventory[customerId][period] + irp.customers[customerId - 1].demand[period]);

        // Visit the nearest customer
        double travelCost = irp.costMatrix[vehicles[nearestVehicleIndex].currentLocation][customerId];
        vehicles[nearestVehicleIndex].visitCustomer(customerId, deliveryAmount);
        updateInventory(period, customerId, currentInventory, deliveryAmount);
        candidateList.erase(candidateList.begin() + nearestCustomerIndex);
        customerServed[customerId] = true; // Mark customer as served

        // Update the current period cost and location
        vehicles[nearestVehicleIndex].currentLocation = customerId;

        // Record the route and update route cost
        solution.vehicleRoutes[period][nearestVehicleIndex].route.push_back({customerId, deliveryAmount});
        solution.vehicleRoutes[period][nearestVehicleIndex].routeCost += travelCost;
    }

    // Ensure all vehicles return to the depot
    for (auto& vehicle : vehicles) {
        double returnCost = irp.costMatrix[vehicle.currentLocation][irp.depots[0].id];
        vehicle.returnToDepot(irp.depots[0].id);
        solution.vehicleRoutes[period][vehicle.id].route.push_back({irp.depots[0].id, 0});
        solution.vehicleRoutes[period][vehicle.id].routeCost += returnCost;

        #ifdef DEBUG
        double recalculatedCost = calculateRouteCost(solution.vehicleRoutes[period][vehicle.id]);
        if (std::abs(recalculatedCost - solution.vehicleRoutes[period][vehicle.id].routeCost) > 1e-6) {
            std::cerr << "Discrepancy in buildRoutes cost calculation for vehicle "
                      << vehicle.id << ": Calculated: " << solution.vehicleRoutes[period][vehicle.id].routeCost
                      << ", Recalculated: " << recalculatedCost << std::endl;
        
        }
        #endif
    }

    solution.currentInventory = currentInventory;
    solution.calculateCosts(); // Calculates the total costs (route and inventory)
    return solution;
}

void Solver::twoOpt(Route& route) { // correto
    int n = route.route.size();
    if (n <= 2) return;

    bool improvement = true;
    while (improvement) {
        improvement = false;

        for (int i = 1; i < n - 2; ++i) {
            for (int j = i + 1; j < n - 1; ++j) {
                if (j - i == 1) continue;

                double oldCost = irp.costMatrix[route.route[i - 1].first][route.route[i].first] +
                                 irp.costMatrix[route.route[j].first][route.route[j + 1].first];

                double newCost = irp.costMatrix[route.route[i - 1].first][route.route[j].first] +
                                 irp.costMatrix[route.route[i].first][route.route[j + 1].first];

                double costImprovement = oldCost - newCost;

                if (costImprovement > 0) {
                    #ifdef DEBUG
                    std::cout << "2-Opt - Route before modification:" << std::endl;
                    route.printRoute();
                    std::cout << "2-Opt - Route cost before modification: " << route.routeCost << std::endl;
                    #endif

                    std::reverse(route.route.begin() + i, route.route.begin() + j + 1);
                    route.routeCost -= costImprovement;
                    improvement = true;

                    #ifdef DEBUG
                    std::cout << "2-Opt - Route after modification:" << std::endl;
                    route.printRoute();
                    std::cout << "2-Opt - Route cost after modification: " << route.routeCost << std::endl;
                    double recalculatedCost = calculateRouteCost(route);
                    if (std::abs(recalculatedCost - route.routeCost) > 1e-6) {
                        std::cerr << "Discrepancy in 2-opt cost calculation: Calculated: "
                                  << route.routeCost << ", Recalculated: " << recalculatedCost << std::endl;
                    return;
                    }
                    #endif
                }
            }
        }
    }
}







void Solver::exchangeOneOne(Route& route) { //exchange 11 implementar 21 22
    int n = route.route.size();
    if (n <= 2) return;

    bool improvement = true;
    while (improvement) {
        improvement = false;

        for (int i = 1; i < n - 1; ++i) {
            for (int j = i + 1; j < n - 1; ++j) {
                
                double oldCost1 = irp.costMatrix[route.route[i - 1].first][route.route[i].first] +
                                  irp.costMatrix[route.route[i].first][route.route[i + 1].first];
                double oldCost2 = irp.costMatrix[route.route[j - 1].first][route.route[j].first] +
                                  irp.costMatrix[route.route[j].first][route.route[j + 1].first];

                if (j == i + 1) {
                    oldCost1 = irp.costMatrix[route.route[i - 1].first][route.route[i].first] +
                               irp.costMatrix[route.route[i].first][route.route[j].first] +
                               irp.costMatrix[route.route[j].first][route.route[j + 1].first];
                    oldCost2 = 0;
                }

                double newCost1 = irp.costMatrix[route.route[i - 1].first][route.route[j].first] +
                                  irp.costMatrix[route.route[j].first][route.route[i + 1].first];
                double newCost2 = irp.costMatrix[route.route[j - 1].first][route.route[i].first] +
                                  irp.costMatrix[route.route[i].first][route.route[j + 1].first];

                if (j == i + 1) {
                    newCost1 = irp.costMatrix[route.route[i - 1].first][route.route[j].first] +
                               irp.costMatrix[route.route[j].first][route.route[i].first] +
                               irp.costMatrix[route.route[i].first][route.route[j + 1].first];
                    newCost2 = 0;
                }

                double costImprovement = (oldCost1 + oldCost2) - (newCost1 + newCost2);

                if (costImprovement > 0) {
                    #ifdef DEBUG
                    std::cout << "exchangeOneOne - Route before modification:" << std::endl;
                    route.printRoute();
                    std::cout << "Route Cost: " << route.routeCost << std::endl;
                    #endif

                    std::swap(route.route[i], route.route[j]);
                    route.routeCost -= costImprovement;
                    improvement = true;

                    #ifdef DEBUG
                    std::cout << "exchangeOneOne - Route after modification:" << std::endl;
                    route.printRoute();
                    std::cout << "Route Cost: " << route.routeCost << std::endl;
                    double recalculatedCost = calculateRouteCost(route);
                    if (std::abs(recalculatedCost - route.routeCost) > 1e-6) {
                        std::cerr << "Discrepancy in exchangeOneOne cost calculation: Calculated: "
                                  << route.routeCost << ", Recalculated: " << recalculatedCost << std::endl;
                    exit(0);
                    }
                    #endif
                    break;
                }
            }

            if (improvement) break;
        }
    }
}

void Solver::exchangeTwoTwo(Route& route) { // ok
    int n = route.route.size();
    if (n <= 3) return; // Ensure there are at least two pairs of customers to swap

    bool improvement = true;
    while (improvement) {
        improvement = false;

        for (int i = 1; i < n - 2; ++i) { 
            for (int j = i + 2; j < n - 2; ++j) { // Avoid overlapping pairs
                // Define two pairs of consecutive customers
                auto customer1 = route.route[i];
                auto customer2 = route.route[i + 1];
                auto customer3 = route.route[j];
                auto customer4 = route.route[j + 1];

                // Calculate current cost for both pairs
                double oldCost1 = irp.costMatrix[route.route[i - 1].first][customer1.first] +
                                  irp.costMatrix[customer2.first][route.route[i + 2].first];
                double oldCost2 = irp.costMatrix[route.route[j - 1].first][customer3.first] +
                                  irp.costMatrix[customer4.first][route.route[j + 2].first];

                // Calculate the cost if swapped
                double newCost1 = irp.costMatrix[route.route[i - 1].first][customer3.first] +
                                  irp.costMatrix[customer4.first][route.route[i + 2].first];
                double newCost2 = irp.costMatrix[route.route[j - 1].first][customer1.first] +
                                  irp.costMatrix[customer2.first][route.route[j + 2].first];

                double costImprovement = (oldCost1 + oldCost2) - (newCost1 + newCost2);

                if (costImprovement > 0) {
                    #ifdef DEBUG
                    std::cout << "exchangeTwoTwo - Route before modification:" << std::endl;
                    route.printRoute();
                    std::cout << "Route Cost: " << route.routeCost << std::endl;
                    #endif

                    // Swap pairs
                    std::swap(route.route[i], route.route[j]);
                    std::swap(route.route[i + 1], route.route[j + 1]);
                    route.routeCost -= costImprovement;
                    improvement = true;

                    #ifdef DEBUG
                    std::cout << "exchangeTwoTwo - Route after modification:" << std::endl;
                    route.printRoute();
                    std::cout << "Route Cost: " << route.routeCost << std::endl;
                    double recalculatedCost = calculateRouteCost(route);
                    if (std::abs(recalculatedCost - route.routeCost) > 1e-6) {
                        std::cerr << "Discrepancy in exchangeTwoTwo cost calculation: Calculated: "
                                  << route.routeCost << ", Recalculated: " << recalculatedCost << std::endl;
                    exit(0);
                    }
                    #endif
                    break;
                }
            }
            if (improvement) break;
        }
    }
}

void Solver::exchangeTwoOne(Route& route) { 
    int n = route.route.size();
    if (n <= 3) return; // Ensure there are enough customers for a 2-1 exchange

    bool improvement = true;
    while (improvement) {
        improvement = false;

        for (int i = 1; i < n - 2; ++i) { // Loop for the first pair
            for (int j = 1; j < n - 1; ++j) { // Loop for the single customer
                if (i == j || i + 1 == j) continue; // Skip overlapping elements

                // Define customers in the swap: two consecutive and one single
                auto customer1 = route.route[i];
                auto customer2 = route.route[i + 1];
                auto customer3 = route.route[j];

                // Calculate current costs
                double oldCost1 = irp.costMatrix[route.route[i - 1].first][customer1.first] +
                                  irp.costMatrix[customer2.first][route.route[i + 2].first];
                double oldCost2 = irp.costMatrix[route.route[j - 1].first][customer3.first] +
                                  irp.costMatrix[customer3.first][route.route[j + 1].first];

                // Calculate new costs if swapped
                double newCost1 = irp.costMatrix[route.route[i - 1].first][customer3.first] +
                                  irp.costMatrix[customer3.first][route.route[i + 2].first];
                double newCost2 = irp.costMatrix[route.route[j - 1].first][customer1.first] +
                                  irp.costMatrix[customer2.first][route.route[j + 1].first];

                double costImprovement = (oldCost1 + oldCost2) - (newCost1 + newCost2);

                if (costImprovement > 0) {
                    #ifdef DEBUG
                    std::cout << "exchangeTwoOne - Route before modification:" << std::endl;
                    route.printRoute();
                    std::cout << "Route Cost: " << route.routeCost << std::endl;
                    #endif

                    // Perform the 2-1 exchange
                    route.route.erase(route.route.begin() + i, route.route.begin() + i + 2);
                    route.route.insert(route.route.begin() + i, customer3);
                    route.route.erase(route.route.begin() + j);
                    route.route.insert(route.route.begin() + j, {customer1, customer2});
                    route.routeCost -= costImprovement;
                    improvement = true;

                    #ifdef DEBUG
                    std::cout << "exchangeTwoOne - Route after modification:" << std::endl;
                    route.printRoute();
                    std::cout << "Route Cost: " << route.routeCost << std::endl;
                    double recalculatedCost = calculateRouteCost(route);
                    if (std::abs(recalculatedCost - route.routeCost) > 1e-6) {
                        std::cerr << "Discrepancy in exchangeTwoOne cost calculation: Calculated: "
                                  << route.routeCost << ", Recalculated: " << recalculatedCost << std::endl;
                    exit(0);
                    }
                    #endif
                    break;
                }
            }
            if (improvement) break;
        }
    }
}



void Solver::reinsertion(Route& route) { // reinsertion
    int n = route.route.size();
    if (n <= 2) return;

    bool improvement = true;
    while (improvement) {
        improvement = false;

        for (int i = 1; i < n - 1; ++i) {  // Skip the depot at the start and end
            for (int j = 1; j < n - 1; ++j) {
                if (i == j || i == j + 1 || j == i + 1) continue;

                // Calculate the cost before removing customer i
                double oldCost = irp.costMatrix[route.route[i - 1].first][route.route[i].first] +
                                 irp.costMatrix[route.route[i].first][route.route[i + 1].first] +
                                 irp.costMatrix[route.route[j - 1].first][route.route[j].first];

                // Calculate the cost after inserting customer i at position j
                double newCost = irp.costMatrix[route.route[i - 1].first][route.route[i + 1].first] +
                                 irp.costMatrix[route.route[j - 1].first][route.route[i].first] +
                                 irp.costMatrix[route.route[i].first][route.route[j].first];

                // Calculate the improvement in cost
                double costImprovement = oldCost - newCost;

                if (costImprovement > 0) {
                    // Print the route before modification
                    #ifdef DEBUG
                    std::cout << "Route before modification (cost improvement found):" << std::endl;
                    route.printRoute();
                    #endif

                    // Remove customer i from its original position
                    auto customer = route.route[i];
                    route.route.erase(route.route.begin() + i);

                    // Adjust j if necessary due to shifting caused by erase
                    if (i < j) {
                        j--; // Since we've removed an element before j
                    }

                    // Insert customer i at the new position j
                    route.route.insert(route.route.begin() + j, customer);

                    // Update the route cost
                    route.routeCost -= costImprovement;
                    improvement = true;

                    #ifdef DEBUG
                    std::cout << "reinsertion improvement made:" << std::endl;
                    std::cout << "  Customer moved from index " << i << " to index " << j << std::endl;
                    std::cout << "  Previous Route Cost: " << (route.routeCost + costImprovement) << std::endl;
                    std::cout << "  Improved Route Cost: " << route.routeCost << std::endl;

                    // Print the route after modification
                    std::cout << "Route after modification:" << std::endl;
                    route.printRoute();

                    // Recalculate and verify route cost
                    double recalculatedCost = calculateRouteCost(route);
                    std::cout << "  Recalculated Route Cost (after improvement): " << recalculatedCost << std::endl;
                    if (std::abs(recalculatedCost - route.routeCost) > 1e-6) {
                        std::cerr << "Discrepancy detected after improvement: Calculated: "
                                  << route.routeCost << ", Recalculated: " << recalculatedCost << std::endl;
                    exit(0);
                    }
                    #endif

                    break; // Stop the loop and restart after a successful relocation
                }
            }

            if (improvement) {
                break; // Restart the search after a relocation
            }
        }
    }
}





void Solver::swapOneOne(Route& route1, Route& route2) { // swap 11
    bool improvement = true;
    while (improvement) {
        improvement = false;

        for (int i = 1; i < route1.route.size() - 1; ++i) {
            for (int j = 1; j < route2.route.size() - 1; ++j) {
                auto& customer1 = route1.route[i];
                auto& customer2 = route2.route[j];

                if (route1.remainingCapacity + customer1.second - customer2.second >= 0 &&
                    route2.remainingCapacity + customer2.second - customer1.second >= 0) {

                    double oldCost1 = irp.costMatrix[route1.route[i - 1].first][customer1.first] +
                                      irp.costMatrix[customer1.first][route1.route[i + 1].first];
                    double oldCost2 = irp.costMatrix[route2.route[j - 1].first][customer2.first] +
                                      irp.costMatrix[customer2.first][route2.route[j + 1].first];

                    double newCost1 = irp.costMatrix[route1.route[i - 1].first][customer2.first] +
                                      irp.costMatrix[customer2.first][route1.route[i + 1].first];
                    double newCost2 = irp.costMatrix[route2.route[j - 1].first][customer1.first] +
                                      irp.costMatrix[customer1.first][route2.route[j + 1].first];

                    double costImprovement = (oldCost1 + oldCost2) - (newCost1 + newCost2);

                    if (costImprovement > 0) {
                        #ifdef DEBUG
                        std::cout << "swapOneOne - Route1 before modification:" << std::endl;
                        route1.printRoute();
                        std::cout << "Route cost before modification: " <<route1.routeCost<< std::endl;
                        std::cout << "swapOneOne Routes - Route2 before modification:" << std::endl;
                        route2.printRoute();
                        std::cout << "Route cost before modification: " <<route2.routeCost<< std::endl;
                        #endif

                        std::swap(route1.route[i], route2.route[j]);
                        route1.routeCost -= oldCost1 - newCost1;
                        route2.routeCost -= oldCost2 - newCost2;
                        route1.remainingCapacity += customer1.second - customer2.second;
                        route2.remainingCapacity += customer2.second - customer1.second;
                        improvement = true;

                        #ifdef DEBUG
                        std::cout << "swapOneOne Routes - Route1 after modification:" << std::endl;
                        route1.printRoute();
                        std::cout << "Route cost after modification: " <<route1.routeCost<< std::endl;
                        std::cout << "swapOneOne Routes - Route2 after modification:" << std::endl;
                        route2.printRoute();
                        std::cout << "Route cost after modification: " <<route2.routeCost<< std::endl;
                        double recalculatedCost1 = calculateRouteCost(route1);
                        double recalculatedCost2 = calculateRouteCost(route2);
                        if (std::abs(recalculatedCost1 - route1.routeCost) > 1e-6 || 
                            std::abs(recalculatedCost2 - route2.routeCost) > 1e-6) {
                            std::cerr << "Discrepancy in swapOneOne cost calculation: "
                                      << "Route1 Calculated: " << route1.routeCost
                                      << ", Recalculated: " << recalculatedCost1 << std::endl;
                            std::cerr << "Route2 Calculated: " << route2.routeCost
                                      << ", Recalculated: " << recalculatedCost2 << std::endl;
                        exit(0);
                        }
                        #endif
                    }
                }
            }
        }
    }
}



void Solver::shiftOneZero(Route& route1, Route& route2) { // shift 10
    bool improvement = true;
   
    // Verificação dos tamanhos para evitar loops desnecessários
    if (route1.route.size() < 2 || route2.route.size() < 2) return;

    while (improvement) {
        improvement = false;

        // Loop sobre os clientes da rota 1 (exceto depósitos nas extremidades)
        for (int i = 1; i < route1.route.size() - 1; ++i) {
            auto customer = route1.route[i];

            // Verifica se a rota 2 tem capacidade suficiente
            for (int j = 1; j < route2.route.size(); ++j) {
                if (route2.remainingCapacity >= customer.second) {
                    // Calcula os custos antigos e novos somente se os índices são válidos
                    if (i - 1 >= 0 && i + 1 < route1.route.size() &&
                        j - 1 >= 0 && j < route2.route.size()) {
                        
                        double oldCost1 = irp.costMatrix[route1.route[i - 1].first][customer.first] +
                                          irp.costMatrix[customer.first][route1.route[i + 1].first];

                        double newCost1 = irp.costMatrix[route1.route[i - 1].first][route1.route[i + 1].first];

                        double oldCost2 = irp.costMatrix[route2.route[j - 1].first][route2.route[j].first];

                        double newCost2 = irp.costMatrix[route2.route[j - 1].first][customer.first] +
                                          irp.costMatrix[customer.first][route2.route[j].first];

                        double costImprovement = (oldCost1 + oldCost2) - (newCost1 + newCost2);

                        if (costImprovement > 0) {
                            #ifdef DEBUG
                            std::cout << "shift 1 0 Routes - Route1 before modification:" << std::endl;
                            route1.printRoute();
                            std::cout << "Route cost before modification: " <<route1.routeCost<< std::endl;
                            std::cout << "shift 1 0 Routes - Route2 before modification:" << std::endl;
                            route2.printRoute();
                            std::cout << "Route cost before modification: " <<route2.routeCost<< std::endl;
                            #endif
                            // Remove o cliente da rota 1 e insere na posição correta da rota 2
                            route1.route.erase(route1.route.begin() + i);
                            route2.route.insert(route2.route.begin() + j, customer);
                            route1.routeCost -= oldCost1 - newCost1;
                            route2.routeCost -= oldCost2 - newCost2;
                            route1.remainingCapacity += customer.second;
                            route2.remainingCapacity -= customer.second;
                            improvement = true;
                            
                        #ifdef DEBUG
                        std::cout << "shift 1 0 Routes - Route1 after modification:" << std::endl;
                        route1.printRoute();
                        std::cout << "Route cost after modification: " <<route1.routeCost<< std::endl;
                        std::cout << "shift 1 0 Routes - Route2 after modification:" << std::endl;
                        route2.printRoute();
                        std::cout << "Route cost after modification: " <<route2.routeCost<< std::endl;

                        double recalculatedCost1 = calculateRouteCost(route1);
                        double recalculatedCost2 = calculateRouteCost(route2);
                        if (std::abs(recalculatedCost1 - route1.routeCost) > 1e-6 || 
                            std::abs(recalculatedCost2 - route2.routeCost) > 1e-6) {
                            std::cerr << "Discrepancy in shiftOneZero cost calculation: "
                                      << "Route 1 Calculated: " << route1.routeCost
                                      << ", Recalculated: " << recalculatedCost1 << std::endl;
                            std::cerr << "Route 2 Calculated: " << route2.routeCost
                                      << ", Recalculated: " << recalculatedCost2 << std::endl;
                        return;
                        }
                        #endif
                            
                            
                            
                            break;
                        }
                    }
                }
            }
            if (improvement) break; // reinicia a busca após relocação bem-sucedida
        }
    }

}


void Solver::swapRoutes(Route& route1, Route& route2) { // excluir
    bool improvement = true;
    while (improvement) {
        improvement = false;

        for (int i = 1; i < route1.route.size() - 1; ++i) {
            for (int j = 1; j < route2.route.size() - 1; ++j) {
                auto& customer1 = route1.route[i];
                auto& customer2 = route2.route[j];

                if (route1.remainingCapacity + customer1.second - customer2.second >= 0 &&
                    route2.remainingCapacity + customer2.second - customer1.second >= 0) {

                    double oldCost1 = irp.costMatrix[route1.route[i - 1].first][customer1.first] +
                                      irp.costMatrix[customer1.first][route1.route[i + 1].first];

                    double oldCost2 = irp.costMatrix[route2.route[j - 1].first][customer2.first] +
                                      irp.costMatrix[customer2.first][route2.route[j + 1].first];

                    double newCost1 = irp.costMatrix[route1.route[i - 1].first][customer2.first] +
                                      irp.costMatrix[customer2.first][route1.route[i + 1].first];

                    double newCost2 = irp.costMatrix[route2.route[j - 1].first][customer1.first] +
                                      irp.costMatrix[customer1.first][route2.route[j + 1].first];

                    double costImprovement = (oldCost1 + oldCost2) - (newCost1 + newCost2);

                    if (costImprovement > 0) {
                        std::swap(route1.route[i], route2.route[j]);
                        route1.routeCost -= oldCost1 - newCost1;
                        route2.routeCost -= oldCost2 - newCost2;
                        route1.remainingCapacity += customer1.second - customer2.second;
                        route2.remainingCapacity += customer2.second - customer1.second;
                        improvement = true;

                        #ifdef DEBUG
                        double recalculatedCost1 = calculateRouteCost(route1);
                        double recalculatedCost2 = calculateRouteCost(route2);
                        if (std::abs(recalculatedCost1 - route1.routeCost) > 1e-6 || 
                            std::abs(recalculatedCost2 - route2.routeCost) > 1e-6) {
                            std::cerr << "Discrepancy in swapRoutes cost calculation: "
                                      << "Route 1 Calculated: " << route1.routeCost
                                      << ", Recalculated: " << recalculatedCost1 << std::endl;
                            std::cerr << "Route 2 Calculated: " << route2.routeCost
                                      << ", Recalculated: " << recalculatedCost2 << std::endl;
                        return;
                        }
                        #endif
                    }
                }
            }
        }
    }
}

void Solver::shiftTwoZeroInterRoute(Route& route1, Route& route2) {
    bool improvement = true;

    while (improvement) {
        improvement = false;

        for (int i = 1; i < route1.route.size() - 2; ++i) {
            auto customer1 = route1.route[i];
            auto customer2 = route1.route[i + 1];

            // Verifica se route2 tem capacidade para os dois clientes
            if (route2.remainingCapacity >= customer1.second + customer2.second && route1.route.size() > 3) {
                double oldCost1 = irp.costMatrix[route1.route[i - 1].first][customer1.first] +
                                  irp.costMatrix[customer2.first][route1.route[i + 2].first];

                double newCost1 = irp.costMatrix[route1.route[i - 1].first][route1.route[i + 2].first];

                double oldCost2 = irp.costMatrix[route2.route.back().first][customer1.first] +
                                  irp.costMatrix[customer1.first][customer2.first];

                double newCost2 = irp.costMatrix[route2.route.back().first][customer1.first] +
                                  irp.costMatrix[customer1.first][customer2.first];

                double costImprovement = (oldCost1 + oldCost2) - (newCost1 + newCost2);

                if (costImprovement > 0) {
                    #ifdef DEBUGNEW
                    std::cout << "Improvement found in shiftTwoZeroInterRoute. Cost reduced by: " 
                              << costImprovement << std::endl;
                    std::cout << "Route1 before: "; route1.printRoute();
                    std::cout << "Route2 before: "; route2.printRoute();
                    #endif

                    route1.route.erase(route1.route.begin() + i, route1.route.begin() + i + 2);
                    route2.route.insert(route2.route.end() - 1, {customer1, customer2});

                    route1.routeCost -= oldCost1 - newCost1;
                    route2.routeCost -= oldCost2 - newCost2;
                    route1.remainingCapacity += customer1.second + customer2.second;
                    route2.remainingCapacity -= customer1.second + customer2.second;
                    improvement = true;

                    #ifdef DEBUGNEW
                    std::cout << "Route1 after: "; route1.printRoute();
                    std::cout << "Route2 after: "; route2.printRoute();
                    
                    double recalculatedCost1 = calculateRouteCost(route1);
                    double recalculatedCost2 = calculateRouteCost(route2);
                    if (std::abs(recalculatedCost1 - route1.routeCost) > 1e-6 || 
                        std::abs(recalculatedCost2 - route2.routeCost) > 1e-6) {
                        std::cerr << "Discrepancy in shiftTwoZeroInterRoute cost calculation: "
                                  << "Route 1 Calculated: " << route1.routeCost
                                  << ", Recalculated: " << recalculatedCost1 << std::endl;
                        std::cerr << "Route 2 Calculated: " << route2.routeCost
                                  << ", Recalculated: " << recalculatedCost2 << std::endl;
                    return;
                    }
                    #endif
                    break;
                }
            }
        }
    }
}

void Solver::swapTwoOneInterRoute(Route& route1, Route& route2) {
 bool improvement = true;

    while (improvement) {
        improvement = false;

        // Índices baseados no tamanho ATUAL das rotas
        const int n1 = route1.route.size();
        const int n2 = route2.route.size();

        for (int i = 1; i < n1 - 2; ++i) { // Garante que há i+1 e i+2
            for (int j = 1; j < n2 - 1; ++j) { // Garante que há j+1

                // Clientes envolvidos na troca
                auto customer1 = route1.route[i];
                auto customer2 = route1.route[i + 1];
                auto customer3 = route2.route[j];

                // Valida capacidade e tamanho mínimo
                if (route1.remainingCapacity + customer1.second + customer2.second - customer3.second < 0) continue;
                if (route2.remainingCapacity + customer3.second - customer1.second - customer2.second < 0) continue;
                if (route1.route.size() <= 3 || route2.route.size() < 2) continue;

                // Cálculo de custos ANTIGOS (rota original)
                double oldCost1 = 
                    irp.costMatrix[route1.route[i - 1].first][customer1.first] +   // Aresta antes de customer1
                    irp.costMatrix[customer1.first][customer2.first] +              // Aresta entre customer1 e customer2
                    irp.costMatrix[customer2.first][route1.route[i + 2].first];     // Aresta após customer2

                double oldCost2 = 
                    irp.costMatrix[route2.route[j - 1].first][customer3.first] +    // Aresta antes de customer3
                    irp.costMatrix[customer3.first][route2.route[j + 1].first];     // Aresta após customer3

                // Cálculo de custos NOVOS (após troca)
                double newCost1 = 
                    irp.costMatrix[route1.route[i - 1].first][customer3.first] +    // Nova aresta antes de customer3
                    irp.costMatrix[customer3.first][route1.route[i + 2].first];     // Nova aresta após customer3

                double newCost2 = 
                    irp.costMatrix[route2.route[j - 1].first][customer1.first] +    // Nova aresta antes de customer1
                    irp.costMatrix[customer1.first][customer2.first] +              // Aresta entre customer1 e customer2
                    irp.costMatrix[customer2.first][route2.route[j + 1].first];     // Nova aresta após customer2

                double costImprovement = (oldCost1 + oldCost2) - (newCost1 + newCost2);

                if (costImprovement > 0) {
                    #ifdef DEBUGNEW
                    std::cout << "[DEBUG] swapTwoOneInterRoute - Melhoria: " << costImprovement 
                              << " | R1: " << route1.routeCost << " -> " << (route1.routeCost - (oldCost1 - newCost1))
                              << " | R2: " << route2.routeCost << " -> " << (route2.routeCost - (oldCost2 - newCost2)) 
                              << std::endl;
                    #endif

                    // Atualiza as rotas
                    route1.route.erase(route1.route.begin() + i, route1.route.begin() + i + 2);
                    route1.route.insert(route1.route.begin() + i, customer3);
                    route2.route.erase(route2.route.begin() + j);
                    route2.route.insert(route2.route.begin() + j, {customer1, customer2});

                    // Atualiza custos de forma segura
                    route1.routeCost = route1.routeCost - oldCost1 + newCost1; // Forma alternativa para evitar erros
                    route2.routeCost = route2.routeCost - oldCost2 + newCost2;

                    // Atualiza capacidades
                    route1.remainingCapacity += customer1.second + customer2.second - customer3.second;
                    route2.remainingCapacity += customer3.second - customer1.second - customer2.second;

                    improvement = true;

                    #ifdef DEBUGNEW
                    // Verificação imediata
                    double recalculatedCost1 = calculateRouteCost(route1);
                    double recalculatedCost2 = calculateRouteCost(route2);
                    if (std::abs(recalculatedCost1 - route1.routeCost) > 1e-6) {
                        std::cerr << "ERRO: Discrepância na Rota 1!\n"
                                  << "  Esperado: " << recalculatedCost1 
                                  << " | Calculado: " << route1.routeCost 
                                  << " | Diferença: " << (route1.routeCost - recalculatedCost1) << std::endl;
                        exit(1);
                    }
                    #endif

                    break;
                }
            }
            if (improvement) break;
        }
    }
}

void Solver::swapTwoTwoInterRoute(Route& route1, Route& route2) {
 bool improvement = true;

    while (improvement) {
        improvement = false;

        // Tamanhos fixos no início da iteração para evitar loops inválidos
        const int n1 = route1.route.size();
        const int n2 = route2.route.size();

        for (int i = 1; i < n1 - 2; ++i) { // i+2 deve ser válido
            for (int j = 1; j < n2 - 2; ++j) { // j+2 deve ser válido

                // Clientes originais nas posições i, i+1 (rota1) e j, j+1 (rota2)
                const auto& customer1 = route1.route[i];
                const auto& customer2 = route1.route[i + 1];
                const auto& customer3 = route2.route[j];
                const auto& customer4 = route2.route[j + 1];

                // Verifica capacidade
                const int deltaCapacityRoute1 = customer1.second + customer2.second - customer3.second - customer4.second;
                const int deltaCapacityRoute2 = customer3.second + customer4.second - customer1.second - customer2.second;
                if (route1.remainingCapacity + deltaCapacityRoute1 < 0) continue;
                if (route2.remainingCapacity + deltaCapacityRoute2 < 0) continue;

                // Cálculo de custos ANTIGOS
                double oldCost1 = 
                    irp.costMatrix[route1.route[i - 1].first][customer1.first] +   // Aresta antes de i
                    irp.costMatrix[customer1.first][customer2.first] +              // Aresta entre i e i+1
                    irp.costMatrix[customer2.first][route1.route[i + 2].first];     // Aresta após i+1

                double oldCost2 = 
                    irp.costMatrix[route2.route[j - 1].first][customer3.first] +   // Aresta antes de j
                    irp.costMatrix[customer3.first][customer4.first] +              // Aresta entre j e j+1
                    irp.costMatrix[customer4.first][route2.route[j + 2].first];     // Aresta após j+1

                // Cálculo de custos NOVOS
                double newCost1 = 
                    irp.costMatrix[route1.route[i - 1].first][customer3.first] +   // Nova aresta antes de i
                    irp.costMatrix[customer3.first][customer4.first] +              // Nova aresta entre i e i+1
                    irp.costMatrix[customer4.first][route1.route[i + 2].first];     // Nova aresta após i+1

                double newCost2 = 
                    irp.costMatrix[route2.route[j - 1].first][customer1.first] +   // Nova aresta antes de j
                    irp.costMatrix[customer1.first][customer2.first] +              // Nova aresta entre j e j+1
                    irp.costMatrix[customer2.first][route2.route[j + 2].first];     // Nova aresta após j+1

                double costImprovement = (oldCost1 + oldCost2) - (newCost1 + newCost2);

                if (costImprovement > 0) {
                    #ifdef DEBUGNEW
                    std::cout << "[DEBUG] swapTwoTwoInterRoute - Melhoria: " << costImprovement 
                              << " | R1: " << route1.routeCost << " -> " << (route1.routeCost - oldCost1 + newCost1)
                              << " | R2: " << route2.routeCost << " -> " << (route2.routeCost - oldCost2 + newCost2) 
                              << std::endl;
                    #endif

                    // Atualiza as rotas (swap seguro)
                    std::swap(route1.route[i], route2.route[j]);
                    std::swap(route1.route[i + 1], route2.route[j + 1]);

                    // Atualiza custos com a fórmula CORRETA: novo_custo = custo_antigo - oldCost + newCost
                    route1.routeCost = route1.routeCost - oldCost1 + newCost1;
                    route2.routeCost = route2.routeCost - oldCost2 + newCost2;

                    // Atualiza capacidades
                    route1.remainingCapacity += deltaCapacityRoute1;
                    route2.remainingCapacity += deltaCapacityRoute2;

                    improvement = true;

                    #ifdef DEBUGNEW
                    // Verificação imediata
                    double recalculatedCost1 = calculateRouteCost(route1);
                    double recalculatedCost2 = calculateRouteCost(route2);
                    if (std::abs(recalculatedCost1 - route1.routeCost) > 1e-6 || 
                        std::abs(recalculatedCost2 - route2.routeCost) > 1e-6) {
                        std::cerr << "ERRO: Discrepância em swapTwoTwoInterRoute!\n"
                                  << "  Rota 1: Esperado=" << recalculatedCost1 
                                  << " | Calculado=" << route1.routeCost 
                                  << "\n  Rota 2: Esperado=" << recalculatedCost2 
                                  << " | Calculado=" << route2.routeCost << std::endl;
                        exit(1);
                    }
                    #endif

                    break;
                }
            }
            if (improvement) break;
        }
    }
}


Solution Solver::findBestSolution(int n, int maxDmax,unsigned seed, int maxNoImp) {
    std::mt19937 rng(seed);  // Inicializa o gerador de números aleatórios com a semente
    std::uniform_int_distribution<int> dmaxDist(1, irp.nPeriods);  // Distribuição de dmax entre 1 e irp.nPeriods

    Solution bestSolution(irp);
    double bestCost = std::numeric_limits<double>::max();
    int bestDmax = 1;  // Para armazenar o melhor dmax
    float bestAlpha = 0.0;  // Para armazenar o melhor alpha

    // Cria o vetor de alphas variando de 0 até 1, incrementando 0.05
    std::vector<float> alphaValues;
    for (float alpha = 0.0; alpha <= 1.0; alpha += 0.05) {
        alphaValues.push_back(alpha);
    }

    // Itera n vezes
    for (int iteration = 0; iteration < n; ++iteration) {
        // Escolhe aleatoriamente o dmax entre 1 e irp.nPeriods
        int dmax = dmaxDist(rng);

        // Embaralha os valores de alpha para garantir diversidade
        std::shuffle(alphaValues.begin(), alphaValues.end(), rng);

        // Itera por cada valor de alpha
        for (float alpha : alphaValues) {
            // Gera uma solução inicial com os valores atuais de dmax e alpha
            Solution solution = solve(dmax, alpha);
            solution = localSearchRandomized(solution, maxNoImp, seed);

            double currentCost = solution.routeCost + solution.inventoryCost;

            // Atualiza a melhor solução se houver melhoria
            if (currentCost < bestCost) {
                bestCost = currentCost;
                bestSolution = solution;
                bestDmax = dmax;   // Atualiza o melhor dmax
                bestAlpha = alpha; // Atualiza o melhor alpha

                // Mensagem de debug opcional
              //  #ifdef DEBUG
                std::cout << "Melhoria encontrada! Custo: " << currentCost
                          << ", dmax: " << dmax
                          << ", alpha: " << alpha << std::endl;
                //#endif
            }
        }
    }

    // Atribui os valores do melhor dmax e alpha encontrados
    bestSolution.bestAlpha = bestAlpha;
    bestSolution.bestDmax = bestDmax;

    // Mensagem final de debug
    #ifdef DEBUG
    std::cout << "Melhor solução encontrada: Custo Total = " << bestCost
              << ", dmax = " << bestDmax
              << ", alpha = " << bestAlpha << std::endl;
    #endif

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



void Solver::executeLocalSearch(int searchIndex, Solution& solution) {
    bool improvement = true;

    while (improvement) {
        improvement = false;

        for (int t = 0; t < irp.nPeriods; ++t) {
            auto& routes = solution.vehicleRoutes[t];

            // Itera sobre todas as rotas no período
            for (size_t i = 0; i < routes.size(); ++i) {
                if (searchIndex <= 3) {
                    // Busca intra-rota
                    switch (searchIndex) {
                        case 1: twoOpt(routes[i]); break;
                        case 2: exchangeOneOne(routes[i]); break;
                        case 3: reinsertion(routes[i]); break;
                    }
                } else if (searchIndex <= 9) {
                    // Busca inter-rotas
                    for (size_t j = i + 1; j < routes.size(); ++j) {
                        Route& r1 = routes[i];
                        Route& r2 = routes[j];
                        double previousCost = r1.routeCost + r2.routeCost;

                        switch (searchIndex) {
                            case 4: swapOneOne(r1, r2); break;
                            case 5: shiftOneZero(r1, r2); break;
                            case 6: swapRoutes(r1, r2); break;
                            case 7: shiftTwoZeroInterRoute(r1, r2); break;
                            case 8: swapTwoOneInterRoute(r1, r2); break;
                            case 9: swapTwoTwoInterRoute(r1, r2); break;
                        }

                        if (r1.routeCost + r2.routeCost < previousCost) {
                            improvement = true;
                            i = 0; // Reinicia a iteração após uma melhoria
                            j = 0;
                        }
                    }
                } else {
                    // Busca local de inventário
                    improvement |= localSearchInventory(solution);
                }
            }
        }
    }
}
// Função principal de busca local randomizada
Solution Solver::localSearchRandomized(Solution& solution, int maxNoImprovement, unsigned int seed) {
    std::mt19937 gen(seed);
    int noImprovementCount = 0;
    double bestCost = solution.inventoryCost + solution.routeCost;

    // Cria um vetor de índices para as buscas locais
    std::vector<int> searchIndices(12); 
    std::iota(searchIndices.begin(), searchIndices.end(), 1);

    while (noImprovementCount < maxNoImprovement) {
        // Embaralha os índices para a ordem de execução das buscas locais
        std::shuffle(searchIndices.begin(), searchIndices.end(), gen);
        bool improvement = false;

        // Aplica cada busca local na ordem embaralhada
        for (int searchIndex : searchIndices) {
            double previousCost = solution.inventoryCost + solution.routeCost;
            executeLocalSearch(searchIndex, solution);

            // Calcula o custo total e verifica se houve melhoria
            double currentCost = solution.inventoryCost + solution.routeCost;
            if (currentCost < previousCost) {
                improvement = true;
                noImprovementCount = 0;
                bestCost = currentCost;
            }
        }

        // Se nenhuma melhoria foi encontrada, incrementa o contador
        if (!improvement) {
            ++noImprovementCount;
        }
    }

    return solution;
}

bool Solver::supplyInsertion(Solution& solution, int periodFrom, int periodTo, int customerId, int quantity) {
    // Validação básica
    if (periodFrom == periodTo || quantity <= 0) return false;
    if (periodFrom < 0 || periodFrom >= irp.nPeriods || periodTo < 0 || periodTo >= irp.nPeriods) return false;
    
    Customer customer = irp.customers[customerId-1];
    auto& inventory = solution.currentInventory;

    // Verifica se a operação é factível
    if (inventory[customerId][periodFrom] < quantity) return false;
    if (inventory[customerId][periodTo] + quantity > customer.maxLevelInv) return false;

    // Encontra veículo com capacidade
    for (Route& route : solution.vehicleRoutes[periodTo]) {
        if (route.remainingCapacity >= quantity) {
            
            // Cria cópia de segurança
            Solution backup = solution;
            
            // Atualiza inventário
            inventory[customerId][periodFrom] -= quantity;
            inventory[customerId][periodTo] += quantity;

            // Atualiza rota
            double oldCost = calculateRouteCost(route);
            route.route.insert(route.route.end()-1, {customerId, quantity});
            route.remainingCapacity -= quantity;
            route.routeCost = calculateRouteCost(route); // Atualização exata
            
            // Atualiza custos de inventário
            solution.inventoryCost -= quantity * customer.invCost * (periodTo - periodFrom);

            // Validação pós-operacional
            if (!validateInventory(solution) || route.remainingCapacity < 0) {
                solution = backup; // Rollback completo
                return false;
            }
            
            return true;
        }
    }
    return false;
}


bool Solver::supplyRemoval(Solution& solution, int period, int customerId, int quantity) {
    if (quantity <= 0) return false;
    auto& inventory = solution.currentInventory;
    
    // Verifica inventário
    if (inventory[customerId][period] < quantity) return false;

    // Procura em todas as rotas do período
    for (Route& route : solution.vehicleRoutes[period]) {
        for (auto it = route.route.begin(); it != route.route.end(); ++it) {
            if (it->first == customerId && it->second >= quantity) {
                
                Solution backup = solution;
                
                // Atualiza inventário
                inventory[customerId][period] -= quantity;
                route.remainingCapacity += quantity;
                
                // Remove entrega
                if (it->second == quantity) {
                    route.route.erase(it);
                } else {
                    it->second -= quantity;
                }
                
                // Atualiza custos
                route.routeCost = calculateRouteCost(route);
                solution.inventoryCost += quantity * irp.customers[customerId-1].invCost;

                // Validação crítica
                if (!validateInventory(solution) || route.remainingCapacity > irp.Capacity) {
                    solution = backup;
                    return false;
                }
                
                return true;
            }
        }
    }
    return false;
}


bool Solver::shiftDelivery(Solution& solution, int periodFrom, int periodTo, int customerId) {
    
    
    // if (periodFrom == periodTo) return false;
    
    // const Customer& customer = irp.customers[customerId - 1];
    // int maxQuantity = solution.currentInventory[customerId][periodFrom];
    // if (maxQuantity <= 0) return false;

    // // Tenta encontrar quantidade máxima transferível
    // for (int q = maxQuantity; q > 0; --q) {
    //     Solution temp = solution;
        
    //     if (supplyRemoval(temp, periodFrom, customerId, q) &&
    //         supplyInsertion(temp, periodFrom, periodTo, customerId, q)) 
    //     {
    //         if (validateSolution(temp)) {
    //             solution = temp;
    //             return true;
    //         }
    //     }
    // }
    // return false;
    
    auto& inventory = solution.currentInventory;
    const Customer& customer = irp.customers[customerId - 1];

    // Validação dos períodos e cliente
    if (periodFrom < 0 || periodFrom >= irp.nPeriods || periodTo < 0 || periodTo >= irp.nPeriods || periodFrom == periodTo)
        return false;
    if (customerId < 1 || customerId > irp.nCustomers) return false;

    int availableQuantity = inventory[customerId][periodFrom];
    if (availableQuantity <= 0) return false;

    // Tenta remover a entrega do período atual
    if (!supplyRemoval(solution, periodFrom, customerId, availableQuantity)) return false;

    // Tenta inserir a entrega no período de destino
    if (!supplyInsertion(solution, periodFrom, periodTo, customerId, availableQuantity)) {
        // Reverte a remoção se falhar
        supplyInsertion(solution, periodTo, periodFrom, customerId, availableQuantity);
        return false;
    }
     return true;
    //return false;
}




// bool Solver::localSearchInventory(Solution& solution) {
//     bool improvement = false;

//     for (int t = 0; t < irp.nPeriods; ++t) {
//         for (auto& customer : irp.customers) {
//             bool localImprovement = true;

//             while (localImprovement) {
//                 localImprovement = false;
//                 double bestCost = solution.inventoryCost + solution.routeCost;

//                 // Armazena o estado inicial da solução
//                 Solution bestSolution = solution;

//                 // std::cout << "\n### Local Search Inventory - Starting Iteration ###" << std::endl;
//                 // std::cout << "Customer ID: " << customer.id << ", Period: " << t 
//                 //           << ", Current Cost: " << bestCost << std::endl;

//                 // Testa supplyInsertion para períodos futuros
//                 for (int periodTo = t + 1; periodTo < irp.nPeriods; ++periodTo) {
//                     for (int quantity = 1; quantity <= irp.Capacity; ++quantity) {
//                         if (supplyInsertion(solution, t, periodTo, customer.id, quantity)) {
//                             double currentCost = solution.inventoryCost + solution.routeCost;
                            
//                             // std::cout << "Customer ID: " << customer.id << ", Period: " << t 
//                             //     << ", Current Cost: " << bestCost << std::endl;
//                             // std::cout << "supplyInsertion: Moved " << quantity 
//                             //           << " from Period " << t << " to " << periodTo
//                             //           << ", New Cost: " << currentCost << std::endl;

//                             if (currentCost < bestCost) {
//                                 bestCost = currentCost;
//                                 bestSolution = solution;  // Guarda a melhor solução
//                                 localImprovement = true;
//                             }

//                             // Reverte a modificação para testar outras opções
//                             supplyRemoval(solution, periodTo, customer.id, quantity);
//                         }
//                     }
//                 }

//                 // Testa supplyRemoval no mesmo período
//                 for (int quantity = 1; quantity <= irp.Capacity; ++quantity) {
//                     if (supplyRemoval(solution, t, customer.id, quantity)) {
//                         double currentCost = solution.inventoryCost + solution.routeCost;

//                         // std::cout << "supplyRemoval: Removed " << quantity 
//                         //           << " from Period " << t
//                         //           << ", New Cost: " << currentCost << std::endl;

//                         if (currentCost < bestCost) {
//                             bestCost = currentCost;
//                             bestSolution = solution;
//                             localImprovement = true;
//                         }

//                         // Reverte a modificação para testar outras opções
//                         supplyInsertion(solution, t, t, customer.id, quantity);
//                     }
//                 }

//                 // Testa shiftDelivery para outros períodos
//                 for (int periodTo = 0; periodTo < irp.nPeriods; ++periodTo) {
//                     if (periodTo == t) continue;  // Evita mover para o mesmo período

//                     if (shiftDelivery(solution, t, periodTo, customer.id)) {
//                         double currentCost = solution.inventoryCost + solution.routeCost;

//                         // std::cout << "shiftDelivery: Shifted Delivery from Period " 
//                         //           << t << " to " << periodTo
//                         //           << ", New Cost: " << currentCost << std::endl;

//                         if (currentCost < bestCost) {
//                             bestCost = currentCost;
//                             bestSolution = solution;
//                             localImprovement = true;
//                         }

//                         // Reverte a modificação se não foi a melhor solução
//                         shiftDelivery(solution, periodTo, t, customer.id);
//                     }
//                 }

//                 // Aplica a melhor modificação encontrada
//                 if (localImprovement) {
//                     solution = bestSolution;
//                     improvement = true;

//                     std::cout << "Improvement Applied: New Total Cost: " 
//                               << solution.inventoryCost + solution.routeCost << std::endl;
//                 }
//             }
//         }
//     }

//     return improvement;
// }

bool Solver::localSearchInventory(Solution& solution) {
    const int MAX_ATTEMPTS = 50;
    bool improvement = false;
    double bestCost = solution.getTotalCost();
    Solution bestSolution = solution;

    std::random_device rd;
    std::mt19937 gen(rd());
    
    for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
        Solution current = solution;
        
        // Escolhe parâmetros aleatórios
        int customerId = rand() % irp.nCustomers + 1;
        int periodFrom = rand() % irp.nPeriods;
        int periodTo = rand() % irp.nPeriods;
        int quantity = rand() % irp.Capacity; // Quantidades menores são mais seguras

        // Escolhe operação aleatória
        int operation = rand() % 3;
        bool success = false;
        
        switch (operation) {
            case 0:
                success = supplyInsertion(current, periodFrom, periodTo, customerId, quantity);
                break;
            case 1:
                success = supplyRemoval(current, periodFrom, customerId, quantity);
                break;
            case 2:
                success = shiftDelivery(current, periodFrom, periodTo, customerId);
                break;
        }

        if (success && validateSolution(current)) {
            double newCost = current.getTotalCost();
            if (newCost < bestCost - 1e-6) { // Tolerância numérica
                bestCost = newCost;
                bestSolution = current;
                improvement = true;
            }
        }
    }

    if (improvement) {
        solution = bestSolution;
        return true;
    }
    return false;
}


bool Solver::validateSolution(Solution& solution) {
    // Verifica níveis de inventário
    for (int c = 1; c <= irp.nCustomers; ++c) {
        for (int t = 0; t < irp.nPeriods; ++t) {
            if (solution.currentInventory[c][t] < 0 || 
                solution.currentInventory[c][t] > irp.customers[c-1].maxLevelInv) {
                return false;
            }
        }
    }

    // Verifica capacidade das rotas
    for (auto& periodRoutes : solution.vehicleRoutes) {
        for (auto& route : periodRoutes) {
            if (route.remainingCapacity < 0) return false;
            if (route.route.empty()) return false;
        }
    }
    
    return true;
}


bool Solver::validateInventory(Solution& solution) {
    for (int c = 1; c <= irp.nCustomers; c++) {
        for (int t = 0; t < irp.nPeriods; t++) {
            int inv = solution.currentInventory[c][t];
            Customer cust = irp.customers[c-1];
            
            // Verifica limites de inventário
            if (inv < cust.minLevelInv || inv > cust.maxLevelInv) {
                return false;
            }
            
            // Verifica continuidade temporal
            if (t > 0) {
                int prev = solution.currentInventory[c][t-1];
                int delta = prev - cust.demand[t-1];
                if (abs(inv - delta) > 1e-6) { // Permite pequenas diferenças numéricas
                    return false;
                }
            }
        }
    }
    return true;
}