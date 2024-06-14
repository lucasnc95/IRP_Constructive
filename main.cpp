
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <climits>
#include "IRP.h"


int main(int argc, char* argv[]) {
    
    IRP irp;
    irp.readDataFromFile(argv[1]);

    // Testar se os dados foram lidos corretamente
    std::cout << "Número de veículos: " << irp.nVehicles << std::endl;
    std::cout << "Número de depósitos: " << irp.nDepots << std::endl;
    std::cout << "Número de clientes: " << irp.nCustomers << std::endl;
    std::cout << "Número de períodos: " << irp.nPeriods << std::endl;

    // Exibir dados do depósito
    std::cout << "Dados do depósito:" << std::endl;
    for (const auto& depot : irp.depots) {
        std::cout << "ID: " << depot.id << ", X: " << depot.x << ", Y: " << depot.y << std::endl;
    }

    // Exibir dados dos clientes
    std::cout << "Dados dos clientes:" << std::endl;
    for (const auto& customer : irp.customers) {
        std::cout << "ID: " << customer.id << ", X: " << customer.x << ", Y: " << customer.y << std::endl;
        for(int i = 0; i < irp.nPeriods ; i++) {
        std::cout<<"Demand "<<i<<" = "<<customer.demand.at(i)<<" ";
        }
        std::cout<<std::endl;
    }

    // Exibir matriz de custos
    std::cout << "Matriz de custos:" << std::endl;
    for (const auto& row : irp.costMatrix) {
        for (const auto& cost : row) {
            std::cout << cost << "\t";
        }
        std::cout << std::endl;
    }

    irp.forwardDeliveryHeuristic();
    irp.printSolution();

    return 0;
}