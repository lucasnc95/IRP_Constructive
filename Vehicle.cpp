#include "Vehicle.h"

Vehicle::Vehicle(int id, int capacity, int startLocation) 
    : id(id), currentCapacity(capacity), currentLocation(startLocation) {
    route.push_back(startLocation); // Start at the depot
}

void Vehicle::visitCustomer(int customerId, int demand) {
    currentCapacity -= demand;
    currentLocation = customerId;
    route.push_back(customerId);
    deliveries.push_back(demand); // Store the delivered quantity
}

void Vehicle::returnToDepot(int depotId) {
    route.push_back(depotId);
    currentLocation = depotId;
}
