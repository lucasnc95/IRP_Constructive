#ifndef IRP_H
#define IRP_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <climits>
#include <random>

class Customer {
public:
    int id;
    double x, y;
    double invCost;
    int initialInv, minLevelInv, maxLevelInv;
    std::vector<int> demand;
};

class Depot {
public:
    int id;
    double x, y;
    double invCost;
    int initialInv, minLevelInv, maxLevelInv;
    std::vector<int> production;
};

class IRP {
public:
    int nVehicles, nDepots, nCustomers, nPeriods;
    int Vehicle_Type, Fleet_Size, Capacity;
    std::vector<Depot> depots;
    std::vector<Customer> customers;
    std::vector<std::vector<int>> costMatrix; // 2D cost matrix

    void readDataFromFile(const std::string& filename);
};

#endif // IRP_H
