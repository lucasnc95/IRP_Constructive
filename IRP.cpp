#include "IRP.h"

void IRP::readDataFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::string aux;
    file >> aux >> nVehicles >> aux >> nDepots >> aux >> nCustomers >> aux >> nPeriods;

    getline(file, line); // Skip line
    getline(file, line); // Skip line
    file >> Vehicle_Type >> Fleet_Size >> Capacity;
    getline(file, line); // Skip line
    getline(file, line); // Skip line

    // Read depot data
    for (int i = 0; i < nDepots; ++i) {
        Depot depot;
        file >> depot.id >> depot.x >> depot.y >> depot.invCost >> depot.initialInv >> depot.minLevelInv >> depot.maxLevelInv;
        depot.production.resize(nPeriods);
        for (int j = 0; j < nPeriods; ++j)
            file >> depot.production[j];
        depots.push_back(depot);
    }

    getline(file, line); // Skip line
    getline(file, line); // Skip line

    // Read customer data
    for (int i = 0; i < nCustomers; ++i) {
        Customer customer;
        file >> customer.id >> customer.x >> customer.y >> customer.invCost >> customer.initialInv >> customer.minLevelInv >> customer.maxLevelInv;
        customer.demand.resize(nPeriods);
        for (int j = 0; j < nPeriods; ++j) {
            file >> customer.demand[j];
        }
        customers.push_back(customer);
    }

    getline(file, line); // Skip line
    getline(file, line); // Skip line
    getline(file, line); // Skip line
    getline(file, line); // Skip line

    // Read cost matrix
    int matrixSize = nDepots + nCustomers;
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
