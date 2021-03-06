#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

#ifndef NIL
#define NIL - (15 << 25)
#endif

#ifndef EPSILON
#define EPSILON 1.0e-3
#endif

using namespace std;

typedef struct {
    unsigned int u; /* edge origin */
    unsigned int v; /* edge destination */
    double w; /* edge weight */
} Edge;

typedef struct {
    Edge e;
    Edge f;
} ConflictingPair;

struct Sets {
    vector <int> p;
    vector <int> rank;
};

bool edgeComparator (Edge e, Edge f) {
    if (e.u == f.u) {
        return e.v < f.v;
    }

    return e.u < f.u;
}

bool edgeWeightComparator (Edge a, Edge b) {
    return (a.w < b.w);
}

bool conflictingPairComparator (ConflictingPair a, ConflictingPair b) {
    if (a.e.u == b.e.u) {
        if (a.e.v == b.e.v) {
            if (a.f.u == b.f.u) {
                return a.f.v < b.f.v;
            }

            return a.f.u < b.f.u;
        }

        return a.e.v < b.e.v;
    }

    return a.e.u < b.e.u;
}

bool areEdgesExtremesEquals (Edge e, Edge f) {
    return (minmax(e.u, e.v) == minmax(f.u, f.v));
}

void SETSinit (Sets * S, unsigned int n) {
    S->p = vector <int> (n, NIL);
    S->rank = vector <int> (n, NIL);

    for (unsigned int i = 0; i < n; i++) {
        S->p[i] = NIL;
        S->rank[i] = NIL;
    }
}

void SETSdestroy (Sets * S) {
    for (unsigned int i = 0; i < S->rank.size(); i++) {
        S->rank[i] = NIL;
    }

    for (unsigned int i = 0; i < S->p.size(); i++) {
        S->p[i] = NIL;
    }
}

void SETSmake (Sets * S, unsigned int x) {
    S->p[x] = x;
    S->rank[x] = 0;
}

int SETSfind (Sets * S, unsigned int x) {
    if (S->p[x] < 0 || x != (unsigned int) S->p[x]) {
        S->p[x] = SETSfind (S, S->p[x]);
    }

    return S->p[x];
}

void SETSlink (Sets * S, unsigned int x, unsigned int y) {
    if (x != y) {
        if (S->rank[x] > S->rank[y]) {
            S->p[y] = x;
        } else {
            S->p[x] = y;

            if (S->rank[x] == S->rank[y]) {
                S->rank[y]++;
            }
        }
    }
}

void SETSunion (Sets * S, unsigned int x, unsigned int y) {
    SETSlink (S, SETSfind (S, x), SETSfind (S, y));
}

Sets connectedComponents (unsigned int n, set <Edge, bool (*) (Edge, Edge)> E) {
    Sets S;

    SETSinit(&S, n);

    /* For each vertex v ∈ V */
    for (unsigned int v = 0; v < n; v++) {
        SETSmake(&S, v);
    }

    /* For each edge (u, v) ∈ E */
    for (set <Edge>::iterator it = E.begin(); it != E.end(); it++) {
        Edge e = (*it);

        if (SETSfind(&S, e.u) != SETSfind(&S, e.v)) {
            SETSunion(&S, e.u, e.v);
        }
    }

    return S;
}

bool sameComponent (Sets S, unsigned int u, unsigned int v) {
    if (SETSfind(&S, u) == SETSfind(&S, v)) {
        return true;
    }

    return false;
}

bool isConnected (unsigned int n, set <Edge, bool (*) (Edge, Edge)> E) {
    bool result = true;

    Sets components = connectedComponents(n, E);

    for (unsigned int u = 0; u < n - 1 && result; u++) {
        for (unsigned int v = u + 1; v < n && result; v++) {
            if (!sameComponent(components, u, v)) {
                result = false;
            }
        }
    }

    SETSdestroy(&components);

    return result;
}

bool readInput (unsigned int * n, set <Edge, bool (*) (Edge, Edge)> * E, 
        set <ConflictingPair,bool (*) (ConflictingPair, ConflictingPair)> * S, 
        char * inputFilePath) {
    ifstream inputFileStream(inputFilePath, ifstream::in);

    if (!inputFileStream.is_open()) {
        return false;
    }

    string line;

    /* Linhas de comentário, iniciadas com o caractere ‘#’. */
    do {
        getline(inputFileStream, line);
    } while (line[0] == '#');

    /* Uma linha com o nome do arquivo. */

    unsigned int m, p;

    /* Uma linha com o número de vértices do grafo (n). */
    inputFileStream >> (*n);

    /* Uma linha com o número de arestas do grafo (m). */
    inputFileStream >> m;

    /* Uma linha com a quantidade de pares de arestas conflitantes (p). */
    inputFileStream >> p;

    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    (*E) = set <Edge, bool (*) (Edge, Edge)> (edgeComparatorPointer);

    bool (*conflictingPairComparatorPointer) (ConflictingPair, ConflictingPair) = 
        conflictingPairComparator;
    (*S) = set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> (
            conflictingPairComparatorPointer);

    /* m linhas, cada uma identificando uma aresta e = (u, v), de custo ce, sendo 0 ≤ u, */
    /* v ≤ n − 1 e 0 ≤ ce ≤ 500, no seguinte formato: u v ce */
    for (unsigned int j = 0; j < m; j++) {
        int u, v;
        double w;

        inputFileStream >> u >> v >> w;

        Edge e;

        if (u < v) {
            e.u = u;
            e.v = v;
        } else {
            e.u = v;
            e.v = u;
        }

        e.w = w;

        (*E).insert(e);
    }

    /* p linhas, cada uma identificando um par de arestas conflitantes {e, f}, sendo */
    /* e = (u, v), f = (x, y) e 0 ≤ u, v, x, y ≤ n − 1, no seguinte formato: u v x y */
    for (unsigned int j = 0; j < p; j++) {
        Edge e, f;
        unsigned int u, v;

        inputFileStream >> u >> v;

        if (u < v) {
            e.u = u;
            e.v = v;
        } else {
            e.u = v;
            e.v = u;
        }

        inputFileStream >> u >> v;

        if (u < v) {
            f.u = u;
            f.v = v;
        } else {
            f.u = v;
            f.v = u;
        }

        set <Edge, bool (*) (Edge, Edge)> :: iterator it = (*E).find(e);

        if (it != (*E).end()) {
            e.w = it->w;
        }

        it = (*E).find(f);

        if (it != (*E).end()) {
            f.w = it->w;
        }

        ConflictingPair cp;

        if (edgeComparator(e, f)) {
            cp.e = e;
            cp.f = f;
        } else {
            cp.e = f;
            cp.f = e;
        }

        (*S).insert(cp);
    }

    return true;
}

bool isFeasible (unsigned int n, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        set <Edge, bool (*) (Edge, Edge)> solution) {
    if (solution.size() != n - 1) {
        return false;
    }

    if (!isConnected(n, solution)) {
        return false;
    }

    for (set <ConflictingPair>::iterator it = S.begin(); it != S.end(); it++) {
        unsigned int counter = 0;

        set <Edge>::iterator it2 = solution.find(it->e);

        if (it2 != solution.end()) {
            counter++;
        }

        it2 = solution.find(it->f);

        if (it2 != solution.end()) {
            counter++;
        }

        if (counter >= 2) {
            return false;
        }
    }

    return true;
}

bool readResult (char * resultFilePath, unsigned int n, double * bestDualBoundValue, 
        int * bestDualBoundIteration, int * totalIterations, double * bestPrimalBoundValue, 
        int * bestPrimalBoundIteration, set <Edge, bool (*) (Edge, Edge)> * bestPrimalSolution) {
    ifstream resultFileStream(resultFilePath, ifstream::in);

    if (!resultFileStream.is_open()) {
        return false;
    }

    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    (*bestPrimalSolution) = set <Edge, bool (*) (Edge, Edge)> (edgeComparatorPointer);

    /* Uma linha com o melhor limitante dual encontrado impresso com 6 casas decimais. */
    resultFileStream >> (*bestDualBoundValue);

    /* Uma linha com a iteração do método do subgradiente (MS) */
    /* em que o melhor limitante dual */
    /* foi encontrado – considere que as iterações do MS são numeradas a partir de 0. */
    resultFileStream >> (*bestDualBoundIteration);

    /* Uma linha com a quantidade de iterações do MS executadas. */
    resultFileStream >> (*totalIterations);

    /* Uma linha com o melhor limitante primal encontrado. */
    resultFileStream >> (*bestPrimalBoundValue);

    /* Uma linha com a iteração do MS em que o melhor */
    /* limitante primal foi encontrado – se o */
    /* melhor limitante primal encontrado foi calculado antes da execução do MS, */
    /* o número −1 deve ser impresso nesta linha. */
    resultFileStream >> (*bestPrimalBoundIteration);

    /* n−1 linhas, sendo n o número de vértices do grafo da instância de teste, */
    /* representando a árvore geradora correspondente ao melhor limitante primal encontrado, */
    /* cada linha identificando uma aresta e = (u, v) presente na árvore, */
    /* sendo 0 ≤ u, v ≤ n − 1, no seguinte formato: u v */
    for (unsigned int i = 0; i < n - 1; i++) {
        Edge e;

        resultFileStream >> e.u >> e.v;

        (*bestPrimalSolution).insert(e);
    }

    return true;
}

bool isResultValid (unsigned int n, set <Edge, bool (*) (Edge, Edge)> E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        double bestDualBoundValue, int bestDualBoundIteration, int totalIterations, 
        double bestPrimalBoundValue, int bestPrimalBoundIteration, 
        set <Edge, bool (*) (Edge, Edge)> bestPrimalSolution) {
    if (bestDualBoundIteration > totalIterations) {
        return false;
    }

    if (bestPrimalBoundIteration > totalIterations) {
        return false;
    }

    if (bestPrimalBoundValue < bestDualBoundValue) {
        return false;
    }

    if (!isFeasible(n, S, bestPrimalSolution)) {
        return false;
    }

    double bestPrimalSolutionValue = 0.0;

    for (set <Edge>::iterator it = bestPrimalSolution.begin(); it != bestPrimalSolution.end(); 
            it++) {
        Edge e = (*it);

        set <Edge>::iterator it2 = E.find(e);

        if (it2 != E.end()) {
            Edge f = (*it2);

            bestPrimalSolutionValue += f.w;
        }
    }

    if (fabs(bestPrimalBoundValue - bestPrimalSolutionValue) > EPSILON) {
        return false;
    }

    return true;
}

bool writeValidation (char * validationFilePath, bool isValid) {
    ofstream validationFileStream(validationFilePath, ofstream::out);

    if (!validationFileStream.is_open()) {
        return false;
    }

    if (isValid) {
        validationFileStream << 0;
    } else {
        validationFileStream << 1;
    }

    validationFileStream << endl;

    return true;
}

int main (int argc, char * argv[]) {
    if (argc != 4) {
        cerr << "Invalid arguments!" << endl;
        return 1;
    }

    unsigned int n;
    set <Edge, bool (*) (Edge, Edge)> E;
    set <ConflictingPair,bool (*) (ConflictingPair, ConflictingPair)> S;

    if (!readInput(&n, &E, &S, argv[1])) {
        cerr << "Error while reading input!" << endl;
    }

    double bestDualBoundValue, bestPrimalBoundValue;
    int bestDualBoundIteration, totalIterations, bestPrimalBoundIteration;
    set <Edge, bool (*) (Edge, Edge)> bestPrimalSolution;

    if (!readResult(argv[2], n, &bestDualBoundValue, &bestDualBoundIteration, &totalIterations, 
                &bestPrimalBoundValue, &bestPrimalBoundIteration, &bestPrimalSolution)) {
        cerr << "Error while reading result!" << endl;
    }

    if (!writeValidation(argv[3], isResultValid(n, E, S, bestDualBoundValue, 
                    bestDualBoundIteration, totalIterations, bestPrimalBoundValue, 
                    bestPrimalBoundIteration, bestPrimalSolution))) {
        cerr << "Error while writing validation!" << endl;
    }

    return 0;
}

