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

    int dmax = 6;//irp.nPeriods; // Define dmax as needed
    int n = 5000; // Number of solutions to generate
    Solution bestSolution = solver.findBestSolution(n, dmax);
    bestSolution.printSolution();
   // bestSolution.calculateCosts();
    double totalRouteCost = bestSolution.routeCost;
    double totalInventoryCost = bestSolution.inventoryCost;
    double totalCost = totalRouteCost + totalInventoryCost;
    std::cout << "Custo Total das Rotas: " << totalRouteCost << std::endl;
    std::cout << "Custo Total do Inventário: " << totalInventoryCost << std::endl;
    std::cout << "Custo Total: " << totalCost << std::endl;
    bool factivel = bestSolution.isFeasible();
    std::cout << "A solução encontrada é factível? " << (factivel ? "Sim" : "Não") << std::endl;
   // bestSolution.printInventoryLevels(irp.nPeriods-1);
    std::string outputFilename = "solution.txt";
    solver.saveSolutionToFile(bestSolution, outputFilename);

    return 0;
}
