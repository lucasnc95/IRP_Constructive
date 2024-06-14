#ifndef IRP_H
#define IRP_H
#include <iostream>
#include <vector>
#include "Location.h"

class IRP{
public:
    int nVehicles;
    int nDepots;
    int nCustomers;
    int nPeriods;
    int Vehicle_Type;
    int Fleet_Size;
    int Capacity;
    std::vector<Depot> depots;
    std::vector<Customer> customers;
    std::vector<std::vector<int>> costMatrix;
    std::vector<std::vector<std::vector<int>>> solution; // Para armazenar as rotas para cada período

    // Função para ler os dados do arquivo .txt
    void readDataFromFile(const std::string& filename);
    
    void updateInventory(int period, Customer& customer);

    void addCustomerToRoute(std::vector<int>& route, int customerId);

    int calculateInsertionCost1(const std::vector<int>& route, int customerId) ;

    int calculateInsertionCost2(const std::vector<int>& route, int customerId, double gamma) ;

    void buildRoutes(int period, std::vector<int>& route, std::vector<int>& candidateList, int& availableVehicles, int& remainingCapacity) ;

    void forwardDeliveryHeuristic() ;

    void init_solution();

    void printSolution();

    void printInventoryLevels(int period);
};


#endif