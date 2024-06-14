#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <climits>
#include <vector>
#include <limits>
#include <cmath>
#include <random>
#include "IRP.h"

    // Função para ler os dados do arquivo .txt
    void IRP::readDataFromFile(const std::string& filename) {
        std::ifstream file(filename);
        std::string line;
        std::string aux;
        // Ler a primeira linha para obter os números de veículos, depósitos, clientes e períodos
        file >> aux >> nVehicles >> aux >> nDepots >> aux >> nCustomers >> aux >> nPeriods;

        // Ignorar as próximas 2 linhas
        getline(file, line);
        getline(file, line);
        file >> Vehicle_Type >> Fleet_Size >> Capacity;
        getline(file, line);
        getline(file, line);
        // Ler os dados dos depósitos
        for (int i = 0; i < nDepots; ++i) {
            Depot depot;
            file >> depot.id >> depot.x >> depot.y >> depot.invCost >> depot.initialInv >> depot.minLevelInv >> depot.maxLevelInv;
            depot.production.resize(nPeriods);
            for (int j = 0; j < nPeriods; ++j)
                file >> depot.production[j];
            depots.push_back(depot);
        }
            getline(file, line);
            getline(file, line);
        // Ler os dados dos clientes
        for (int i = 0; i < nCustomers; ++i) {
            Customer customer;
            file >> customer.id >> customer.x >> customer.y >> customer.invCost >> customer.initialInv >> customer.minLevelInv >> customer.maxLevelInv;
            customer.demand.resize(nPeriods);
            customer.currentInventory.resize(nPeriods); // Inicializar o vetor de inventário atual
            customer.currentInventory.at(0)= customer.initialInv;
            for (int j = 0; j < nPeriods; ++j) {
                file >> customer.demand[j];
               
            }
           // customer.currentInventory.push_back(customer.initialInv);
            customers.push_back(customer);
        }
            getline(file, line);
            getline(file, line);
            getline(file, line);
            getline(file, line);
        
        // Ler a matriz de custos
        costMatrix.resize(nDepots + nCustomers , std::vector<int>(nDepots + nCustomers , 0));
        int contAux = 0;
        for (int i = 0; i < nDepots + nCustomers ; i++) {
            if (i > 0)
            file >> aux;

            for (int j = 0; j <= contAux ; j++) {
                
                if (j == i)
                    costMatrix[i][j] == 0;
                else{
                    file >> costMatrix[i][j];
                    costMatrix[j][i] = costMatrix[i][j];
                     
                    }
               
            }
            
       contAux++;     
      
        }
    }

    void IRP::init_solution(){

        solution.resize(nPeriods);
        for (int t = 0; t < nPeriods; ++t) {
            solution[t].resize(0);

    //    for (auto customer: customers)
           // customer.currentInventory.push_back(customer.initialInv);
      

    }}


 void IRP::updateInventory(int period, Customer& customer) {
    
    std::cout << "Period " << period << ":\n";
   // std::cout << "CurrentInventory size " << customer.currentInventory.size() << ":\n";
    std::cout<<"Customer "<<customer.id<< " currentInventory value: " << customer.currentInventory.at(period) << ":\n";
    int newInventory = customer.currentInventory.at(period) - customer.demand[period]; 
    customer.currentInventory.at(period) = newInventory;
    std::cout <<"Customer "<<customer.id<< " currentInventory value after demand " << newInventory << ":\n";
   // std::cout << "End period " << period << ":\n";
    
}

    

    void IRP::addCustomerToRoute(std::vector<int>& route, int customerId) {
        route.push_back(customerId);
    }

    int IRP::calculateInsertionCost1(const std::vector<int>& route, int customerId) {
        double minDistance = std::numeric_limits<double>::max();
        int bestPosition = -1;
        for (int i = 0; i < route.size(); ++i) {
            int lastCustomer = route[i];
            double distance = costMatrix[lastCustomer][customerId];
            if (distance < minDistance) {
                minDistance = distance;
                bestPosition = i;
            }
        }
        return bestPosition;
    }

    int IRP::calculateInsertionCost2(const std::vector<int>& route, int customerId, double gamma) {
        double minCost = std::numeric_limits<double>::max();
        int bestPosition = -1;
        for (int i = 0; i < route.size() - 1; ++i) {
            int c1 = route[i];
            int c2 = route[i + 1];
            double cost = costMatrix[c1][customerId] + costMatrix[customerId][c2] - costMatrix[c1][c2] + gamma;
            if (cost < minCost) {
                minCost = cost;
                bestPosition = i + 1;
            }
        }
        return bestPosition;
    }

    void IRP::buildRoutes(int period, std::vector<int>& route, std::vector<int>& candidateList, int& availableVehicles, int& remainingCapacity) {
      std::cout<<"br  "<<std::endl;
       while (!candidateList.empty() && availableVehicles > 0) {
        int customerId = candidateList.back();
        candidateList.pop_back();
        std::cout<<"1  "<<std::endl;
        Customer& customer = customers[customerId];
        int demand = customer.demand.at(period);
        std::cout<<"2 "<<std::endl;
        if (remainingCapacity < demand) {
            --availableVehicles;
            if (availableVehicles > 0) {
                remainingCapacity = Capacity;
            } else {
                break;
            }
        }std::cout<<"3 "<<std::endl;

        if (remainingCapacity >= demand && customer.currentInventory[period] + demand <= customer.maxLevelInv) {
            addCustomerToRoute(route, customerId);
            remainingCapacity -= demand;
            customer.currentInventory[period] += demand;
        }std::cout<<"4 "<<std::endl;
    }std::cout<<"F br  "<<std::endl;
    }

    void IRP::forwardDeliveryHeuristic() {
        init_solution();
        std::cout<<"FD  "<<std::endl;
       for (int t = 0; t < nPeriods; ++t) {
        std::vector<int> candidateList;
        for (auto& customer : customers) {
            updateInventory(t, customer); // Atualiza o estoque para todos os clientes
            if (customer.currentInventory[t] < customer.demand[t]) {
                candidateList.push_back(customer.id);
            }
        }

        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(candidateList.begin(), candidateList.end(), g);

        int availableVehicles = Fleet_Size;
        int remainingCapacity = Capacity;

        while (!candidateList.empty() && availableVehicles > 0) {
            std::vector<int> route = { depots[0].id };
            buildRoutes(t, route, candidateList, availableVehicles, remainingCapacity);
            if (route.size() > 1) { // Verifica se a rota contém mais do que apenas o depósito
                solution[t].push_back(route);
            }
        }
    }
    }

    void IRP::printSolution() {
        for (int t = 0; t < nPeriods; ++t) {
            std::cout << "Period " << t << ":\n";
            for (const auto& route : solution[t]) {
                std::cout << "  Route: ";
                for (int customerId : route) {
                    std::cout << customerId << " ";
                }
                std::cout << "\n";
            }
          //  printInventoryLevels(t);
        }
    }


    void IRP::printInventoryLevels(int period) {
        std::cout << "Inventory levels at period " << period << ":\n";
        for (const auto& depot : depots) {
            std::cout << "  Depot " << depot.id << ": " << depot.currentInventory[period] << "\n";
        }
        for (const auto& customer : customers) {
            std::cout << "  Customer " << customer.id << ": " << customer.currentInventory[period] << "\n";
        }
    }
