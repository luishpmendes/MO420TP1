#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

bool readResult (char * resultFilePath, double * bestDualBoundValue, 
        double * bestPrimalBoundValue) {
    int bestDualBoundIteration, totalIterations;
    ifstream resultFileStream(resultFilePath, ifstream::in);

    if (!resultFileStream.is_open()) {
        return false;
    }

    /* Uma linha com o melhor limitante dual encontrado impresso com 6 casas decimais. */
    resultFileStream >> (*bestDualBoundValue);

    /* Uma linha com a iteração do método do subgradiente (MS) */
    /* em que o melhor limitante dual */
    /* foi encontrado – considere que as iterações do MS são numeradas a partir de 0. */
    resultFileStream >> bestDualBoundIteration;

    /* Uma linha com a quantidade de iterações do MS executadas. */
    resultFileStream >> totalIterations;

    /* Uma linha com o melhor limitante primal encontrado. */
    resultFileStream >> (*bestPrimalBoundValue);

    return true;
}

int main (int argc, char * argv[]) {
    if (argc < 2) {
        cerr << "Invalid arguments!" << endl;
        return 1;
    }

    double meanGap = 0.0;

    for (int i = 1; i < argc; i++) {
        double bestDualBoundValue, bestPrimalBoundValue, gap;
        
        if (!readResult(argv[i], &bestDualBoundValue, &bestPrimalBoundValue)) {
            cerr << "Error while reading result!" << endl;
            return 1;
        }

        gap = (bestPrimalBoundValue - ceil(bestDualBoundValue)) / ceil(bestDualBoundValue);
        meanGap += gap;

        cout << argv[i] << " ";
        cout << ((int) bestPrimalBoundValue) << " ";
        cout << fixed << setprecision(6) << bestDualBoundValue << " ";
        cout << fixed << setprecision(6) << gap << endl;
    }

    meanGap /= ((double) (argc - 1.0));

    cout << fixed << setprecision(6) << meanGap << endl;

    return 0;
}

