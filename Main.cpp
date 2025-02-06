#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include "IRP.h"
#include "Solver.h"

#include <fstream>
#include <string>
#include <iostream>

void saveResultsToCSV(const std::string& instanceName, int dmax, double bestAlpha, 
                      double routeCost, double inventoryCost, double totalCost, 
                      double executionTime, unsigned int seed, const std::string& executableName) {
    std::string filename = "results.csv";
    std::ofstream file;

    // Verifica se o arquivo existe e se está vazio; caso contrário, escreve o cabeçalho
    bool fileExists = std::ifstream(filename).good();
    file.open(filename, std::ios::app);

    if (!fileExists) {
        // Cabeçalho do arquivo CSV
        file << "Nome da Instancia, Melhor dmax, Melhor Alfa, Custo da Rota, Custo do Inventario, "
                "Custo Total, Tempo de Execucao (s), Semente, Nome do Executável\n";
    }

    // Escreve os resultados
    file << instanceName << "," << dmax << "," << bestAlpha << ","
         << routeCost << "," << inventoryCost << "," << totalCost << "," 
         << executionTime << "," << seed << "," << executableName << "\n";

    file.close();
}


int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Uso: " << argv[0] << " <input_file> <n> <seed> <maxNoImprovement>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string instanceName = argv[1];
    int n = std::stoi(argv[2]);
    unsigned seed = std::stoi(argv[3]);
    int maxNoImprovement = std::stoi(argv[4]);
    IRP irp;
    irp.readDataFromFile(instanceName);

    Solver solver(irp);

    int dmax = irp.nPeriods; 
    auto start = std::chrono::steady_clock::now(); // Tempo de início

    Solution bestSolution = solver.findBestSolution(n, dmax, seed, maxNoImprovement);

    auto end = std::chrono::steady_clock::now();   // Tempo de término
    std::chrono::duration<double> elapsedSeconds = end - start; // Duração em segundos

    double totalRouteCost = bestSolution.routeCost;
    double totalInventoryCost = bestSolution.inventoryCost;
    double totalCost = totalRouteCost + totalInventoryCost;

    std::cout << "Custo Total das Rotas: " << totalRouteCost << std::endl;
    std::cout << "Custo Total do Inventário: " << totalInventoryCost << std::endl;
    std::cout << "Custo Total: " << totalCost << std::endl;

    
    std::string outputFilename = "solution.txt";
    solver.saveSolutionToFile(bestSolution, outputFilename);

    // Salva os resultados no arquivo .csv
    double bestAlpha = bestSolution.bestAlpha; 
    int bestD = bestSolution.bestDmax;
    bestSolution.printSolution();
    saveResultsToCSV(instanceName, bestD, bestAlpha, totalRouteCost, totalInventoryCost, totalCost, elapsedSeconds.count(), seed, argv[0]);

    return 0;
}



























//bool factivel = bestSolution.isFeasible();
   // std::cout << "A solução encontrada é factível? " << (factivel ? "Sim" : "Não") << std::endl;
