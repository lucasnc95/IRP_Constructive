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

    int dmax = 3; // Define dmax as needed
    int n = 10; // Number of solutions to generate
    Solution bestSolution = solver.findBestSolution(n, dmax);

    std::string outputFilename = "solution.txt";
    solver.saveSolutionToFile(bestSolution, outputFilename);

    return 0;
}
