#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#ifndef INFINITE
#define INFINITE 15 << 25
#endif

using namespace std;

bool readOutput (char * outputFilePath, unsigned int * totalTime, 
        unsigned int * bestPrimalBoundTime, unsigned int * bestDualBoundTime) {
    unsigned int iteration, time, bestPrimalBoundValue;
    string primalBoundValueStr; 
    double dualBoundValue, bestDualBoundValue;
    ifstream outputFileStream(outputFilePath, ifstream::in);

    (*totalTime) = 0;
    (*bestPrimalBoundTime) = 0;
    (*bestDualBoundTime) = 0;

    bestDualBoundValue = -INFINITE;
    bestPrimalBoundValue = INFINITE;

    if (!outputFileStream.is_open()) {
        return false;
    }

    while (outputFileStream >> iteration >> time >> dualBoundValue >> primalBoundValueStr) {
        string nanStr ("NaN");

        if ((*totalTime) < time) {
            (*totalTime) = time;
        }

        if (bestDualBoundValue < dualBoundValue) {
            bestDualBoundValue = dualBoundValue;
            (*bestDualBoundTime) = time;
        }

        if (primalBoundValueStr.compare(nanStr) != 0) {
            unsigned int primalBoundValue = stoul(primalBoundValueStr);

            if (bestPrimalBoundValue > primalBoundValue) {
                bestPrimalBoundValue = primalBoundValue;
                (*bestPrimalBoundTime) = time;
            }
        }
    }

    return true;
}

int main (int argc, char * argv[]) {
    if (argc < 2) {
        cerr << "Invalid arguments!" << endl;
        return 1;
    }

    for (int i = 1; i < argc; i++) {
        unsigned int totalTime, bestPrimalBoundTime, bestDualBoundTime;

        if (!readOutput(argv[i], &totalTime, &bestPrimalBoundTime, &bestDualBoundTime)) {
            cerr << "Error while reading output!" << endl;
            return 1;
        }

        cout << argv[i] << " ";
        cout << totalTime << " ";
        cout << bestPrimalBoundTime << " ";
        cout << bestDualBoundTime << endl;
    }

    return 0;
}

