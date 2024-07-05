#include <iostream>
#include "IRP.h"
#include "Solver.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return EXIT_FAILURE;
    }

    IRP irp;
    irp.readDataFromFile(argv[1]);

    Solver solver(irp);

    std::cout << "Número de veículos: " << irp.nVehicles << std::endl;
    std::cout << "Número de depósitos: " << irp.nDepots << std::endl;
    std::cout << "Número de clientes: " << irp.nCustomers << std::endl;
    std::cout << "Número de períodos: " << irp.nPeriods << std::endl;

    std::cout << "Dados do depósito:" << std::endl;
   
   
    for (const auto& depot : irp.depots) {
        std::cout << "ID: " << depot.id << ", X: " << depot.x << ", Y: " << depot.y << std::endl;
    }

    std::cout << "Dados dos clientes:" << std::endl;
    for (const auto& customer : irp.customers) {
        std::cout << "ID: " << customer.id << ", X: " << customer.x << ", Y: " << customer.y << std::endl;
        for (int i = 0; i < irp.nPeriods; i++) {
            std::cout << "Demand " << i << " = " << customer.demand.at(i) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Matriz de custos:" << std::endl;
    for (int i = 0; i < irp.costMatrix.size(); ++i) {
        for (int j = 0; j < irp.costMatrix[i].size(); ++j) {
            std::cout << irp.costMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }

    int dmax = 12; // Define dmax as needed
    int n = 150; // Number of solutions to generate
    Solution bestSolution = solver.findBestSolution(n, dmax);

    std::cout << "Melhor solução encontrada:" << std::endl;
    bestSolution.printSolution();

    double totalRouteCost = bestSolution.routeCost;
    double totalInventoryCost = bestSolution.inventoryCost;
    double totalCost = totalRouteCost + totalInventoryCost;
    std::cout << "Custo Total das Rotas: " << totalRouteCost << std::endl;
    std::cout << "Custo Total do Inventário: " << totalInventoryCost << std::endl;
    std::cout << "Custo Total: " << totalCost << std::endl;

    std::cout << "Estoque em cada período:" << std::endl;
    for (int t = 0; t < irp.nPeriods; ++t) {
      //  bestSolution.printInventoryLevels(t);
    }

    bool factivel = bestSolution.isFeasible();
    std::cout << "A solução encontrada é factível? " << (factivel ? "Sim" : "Não") << std::endl;

    return 0;
}
