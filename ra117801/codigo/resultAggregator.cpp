#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#ifndef INFINITE
#define INFINITE 15 << 25
#endif

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
    unsigned int validGapCounter = 0;

    for (int i = 1; i < argc; i++) {
        double bestDualBoundValue, bestPrimalBoundValue, gap;
        
        if (!readResult(argv[i], &bestDualBoundValue, &bestPrimalBoundValue)) {
            cerr << "Error while reading result!" << endl;
            return 1;
        }

        if (ceil(bestDualBoundValue) <= -INFINITE) {
            gap = nan("");
        } else if (bestPrimalBoundValue >= INFINITE) {
            gap = INFINITE;
        } else {
            gap = (bestPrimalBoundValue - ceil(bestDualBoundValue)) / ceil(bestDualBoundValue);
            meanGap += gap;
            validGapCounter++;
        }

        cout << argv[i] << " ";
        cout << ((int) bestPrimalBoundValue) << " ";
        cout << fixed << setprecision(6) << bestDualBoundValue << " ";
        if (isnan(gap)) {
            cout << "NaN" << endl;
        } else if (gap >= INFINITE) {
            cout << "INFINITE" << endl;
        } else {
            cout << fixed << setprecision(6) << gap << endl;
        }
    }

    meanGap /= ((double) validGapCounter);

    cout << validGapCounter << endl;
    cout << fixed << setprecision(6) << meanGap << endl;

    return 0;
}

