#ifndef LOCATION_H
#define LOCATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <climits>


class Customer {
public:
    int id;
    double x, y;
    double invCost;
    int initialInv, minLevelInv, maxLevelInv;
    std::vector<int> demand;
    std::vector<int> currentInventory; // Valor atual do inventário para cada período
};


class Depot {
public:
    int id;
    double x, y;
    double invCost;
    int initialInv, minLevelInv, maxLevelInv;
    std::vector<int> production;
    std::vector<int> currentInventory; // Valor atual do inventário para cada período
};

#endif