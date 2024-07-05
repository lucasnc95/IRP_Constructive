#ifndef VEHICLE_H
#define VEHICLE_H

#include <vector>

class Vehicle {
public:
    int id;
    int currentCapacity; // capacidade residual
    int currentLocation; // to do remover (solver) 
    // acrescentar custo da rota
    std::vector<int> route;
    std::vector<int> deliveries; // Quantities delivered to each customer

    Vehicle(int id, int capacity, int startLocation);

    void visitCustomer(int customerId, int demand);
    void returnToDepot(int depotId);
};

#endif // VEHICLE_H
