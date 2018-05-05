#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#ifndef INFINITE
#define INFINITE 15 << 25
#endif

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

bool comparator (Edge a, Edge b) {
    return (a.w < b.w);
}

double kruskal (vector <Edge> * A, unsigned int n, vector <Edge> E) {
    double result = 0.0;
    (*A) = vector <Edge> (); /* A ← ∅ */
    Sets S;

    SETSinit (&S, n);

    /* For each vertex v ∈ V */
    for (unsigned int v = 0; v < n; v++) {
        SETSmake (&S, v);
    }

    /* Sort the edges of E into nondecreasing order by weight w */
    sort(E.begin(), E.end(), comparator);

    /* For each edge (u, v) ∈ E, taken in nondecreasing order by weight */
    for (vector <Edge>::iterator it = E.begin(); it != E.end(); it++) {
        Edge e = (*it);
        if (SETSfind (&S, e.u) != SETSfind (&S, e.v)) {
            (*A).push_back(e);
            result += e.w;
            SETSunion (&S, e.u, e.v);
        }
    }

    SETSdestroy (&S);

    return result;
}

Sets connectedComponents (unsigned int n, vector <Edge> E) {
    Sets S;

    SETSinit(&S, n);

    /* For each vertex v ∈ V */
    for (unsigned int v = 0; v < n; v++) {
        SETSmake(&S, v);
    }

    /* For each edge (u, v) ∈ E */
    for (vector <Edge>::iterator it = E.begin(); it != E.end(); it++) {
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

bool readParameters (unsigned int k, unsigned int * timeLimit, double * pi, unsigned int * N, 
        double * minPi) {
    ifstream parameterFileStream(k == 1 ? "param1" : "param2", ifstream::in);
    if (!parameterFileStream.is_open()) {
        return false;
    }
    /* Time limit in seconds */
    parameterFileStream >> (*timeLimit);
    /* π */
    parameterFileStream >> (*pi);
    /* N */
    parameterFileStream >> (*N);
    /* minPi */
    parameterFileStream >> (*minPi);
    return true;
}

bool readInput (unsigned int * n, vector <Edge> * E, vector <ConflictingPair> * S, 
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
    (*E) = vector <Edge> (m);
    (*S) = vector <ConflictingPair> (p);
    /* m linhas, cada uma identificando uma aresta e = (u, v), de custo ce, sendo 0 ≤ u, */
    /* v ≤ n − 1 e 0 ≤ ce ≤ 500, no seguinte formato: u v ce */
    for (unsigned int j = 0; j < m; j++) {
        inputFileStream >> (*E)[j].u >> (*E)[j].v >> (*E)[j].w;
    }
    /* p linhas, cada uma identificando um par de arestas conflitantes {e, f}, sendo */
    /* e = (u, v), f = (x, y) e 0 ≤ u, v, x, y ≤ n − 1, no seguinte formato: u v x y */
    for (unsigned int j = 0; j < p; j++) {
        inputFileStream >> (*S)[j].e.u >> (*S)[j].e.v >> (*S)[j].f.u >> (*S)[j].f.v;
        for (vector <Edge>::iterator it = (*E).begin(); it != (*E).end(); it++) {
            if ((*S)[j].e.u == it->u && (*S)[j].e.v == it->v) {
                (*S)[j].e.w = it->w;
            }
            if ((*S)[j].f.u == it->u && (*S)[j].f.v == it->v) {
                (*S)[j].f.w = it->w;
            }
        }
    }
    return true;
}

vector <double> initialLagrangeMultipliers (int m) {
    vector <double> result (m, 1.0);
    return result;
}

bool termination (chrono::high_resolution_clock::time_point tBegin, unsigned int timeLimit, 
        double bestDualBoundValue, double bestPrimalBoundValue, double pi, double minPi) {
    chrono::high_resolution_clock::time_point tCurrent = chrono::high_resolution_clock::now();
    chrono::seconds elapsedTime = chrono::duration_cast <chrono::seconds> (tCurrent - tBegin);
    if ((unsigned int) elapsedTime.count() >= timeLimit) {
        return true;
    }
    /* It is clear that the subgradient optimization procedure can be terminated if we find */
    /* that bestPrimalBoundValue == bestDualBoundValue.  */
    /* In this case the value of the maximum lower bound coincides with the value of a feasible */
    /* solution and so it must be optimal */
    if (bestDualBoundValue == bestPrimalBoundValue) {
        return true;
    }
    /* Terminate the subgradient optimization procedure when π is small */
    if (pi <= minPi) {
        return true;
    }
    return false;
}

bool isFeasible (unsigned int n, vector <ConflictingPair> S, vector <Edge> solution) {
    if (solution.size() != n - 1) {
        return false;
    }
    Sets components = connectedComponents(n, solution);
    for (unsigned int u = 0; u < n; u++) {
        for (unsigned int v = u + 1; v < n; v++) {
            if (!sameComponent(components, u, v)) {
                SETSdestroy(&components);
                return false;
            }
        }
    }
    SETSdestroy(&components);
    for (vector <ConflictingPair>::iterator it = S.begin(); it != S.end(); it++) {
        unsigned int counter = 0;
        for (vector <Edge>::iterator it2 = solution.begin(); it2 != solution.end(); it2++) {
            if ((it->e.u == it2->u && it->e.v == it2->v) || 
                    (it->f.u == it2->u && it->f.v == it2->v)) {
                counter++;
            }
            if (counter >= 2) {
                return false;
            }
        }
    }
    return true;
}

double fixSolution (vector <Edge> * primalSolution, unsigned int n, vector <Edge> E, 
        vector <ConflictingPair> S) {
    double primalBoundValue = 0.0;
    /* sort the edges of E into nondecreasing order by weight w */
    sort(E.begin(), E.end(), comparator);
    /* Fixing the weight of solution's edges */
    for (vector <Edge>::iterator it = (*primalSolution).begin(); it != (*primalSolution).end(); 
            it++) {
        for (vector <Edge>::iterator it2 = E.begin(); it2 != E.end(); it2++) {
            if (it->u == it2->u && it->v == it2->v) {
                it->w = it2->w;
            }
        }
        primalBoundValue += it->w;
    }
    while (!isFeasible(n, S, (*primalSolution))) {
        /* Choosing a conflicting edge to remove */
        Edge edgeToRemove;
        unsigned int solutionIndexToRemove, EIndexToRemove;
        bool flag = false;
        for (vector <ConflictingPair>::iterator it = S.begin(); it != S.end() && !flag; it++) {
            unsigned int counter = 0;
            for (vector <Edge>::iterator it2 = (*primalSolution).begin(); 
                    it2 != (*primalSolution).end() && !flag; it2++) {
                if (counter >= 2) {
                    if (it->e.w > it->f.w) {
                        edgeToRemove = it->e;
                    } else {
                        edgeToRemove = it->f;
                    }
                    flag = true;
                } else if ((it->e.u == it2->u && it->e.v == it2->v) || 
                        (it->f.u == it2->u && it->f.v == it2->v)){
                    counter++;
                }
            }
        }
        for (unsigned int i = 0; i < (*primalSolution).size(); i++) {
            if ((*primalSolution)[i].u == edgeToRemove.u 
                    && (*primalSolution)[i].v == edgeToRemove.v) {
                solutionIndexToRemove = i;
                break;
            }
        }
        for (unsigned int i = 0; i < E.size(); i++) {
            if (E[i].u == edgeToRemove.u && E[i].v == edgeToRemove.v) {
                EIndexToRemove = i;
                break;
            }
        }
        /* Removing conflicting edge from solution */
        (*primalSolution).erase((*primalSolution).begin() + solutionIndexToRemove);
        /* Updating solution value */
        primalBoundValue -= edgeToRemove.w;
        /* Removing conflicting edge from graph */
        E.erase(E.begin() + EIndexToRemove);
        /* Choosing a non-conflicting edge that re-connects the solution */
        Edge edgeToInsert;
        Sets components = connectedComponents(n, (*primalSolution));
        bool foundEdge = false;
        /* for each edge (u, v) ∈ E, taken in nondecreasing order by weight */
        for (vector <Edge>::iterator it = E.begin(); it != E.end() && !foundEdge; it++) {
            edgeToInsert = (*it);
            /* Checking if the edge connects the solution */
            if (!sameComponent(components, edgeToInsert.u, edgeToInsert.v)) {
                /* Checking if the edge conflicts if another edge in solution */
                bool conflict = false;
                for (vector <ConflictingPair>::iterator it = S.begin(); it != S.end() && 
                        !conflict; it++) {
                    if ((it->e.u == edgeToInsert.u && it->e.v == edgeToInsert.v) || 
                            (it->f.u == edgeToInsert.u && it->f.v == edgeToInsert.v)) {
                        for (vector <Edge>::iterator it2 = (*primalSolution).begin(); 
                                it2 != (*primalSolution).end(); it2++) {
                            if ((it->e.u == it2->u && it->e.v == it2->v) || 
                                    (it->f.u == it2->u && it->f.v == it2->v)) {
                                conflict = true;
                            }
                        }
                    }
                }
                foundEdge = !conflict;
            }
        }
        /* Adding edge to solution */
        (*primalSolution).push_back(edgeToInsert);
        /* Updating solution value */
        primalBoundValue += edgeToInsert.w;
    }
    return primalBoundValue;
}

bool relaxLag1 (double * bestDualBoundValue, int * bestDualBoundIteration, int * totalIterations, 
        double * bestPrimalBoundValue, int * bestPrimalBoundIteration, 
        vector <Edge> * bestPrimalSolution, unsigned int n, vector <Edge> E, 
        vector <ConflictingPair> S, chrono::high_resolution_clock::time_point tBegin, 
        unsigned int timeLimit, double pi, unsigned int N, double minPi) {
    (*bestDualBoundValue) = -INFINITE;
    (*bestDualBoundIteration) = 0;
    (*totalIterations) = 0;
    (*bestPrimalBoundValue) = INFINITE;
    (*bestPrimalBoundIteration) = -1;
    vector <double> u = initialLagrangeMultipliers(S.size());
    unsigned int iterationsWithoutImprovment = 0;
    while (!termination(tBegin, timeLimit, (*bestDualBoundValue), (*bestPrimalBoundValue), pi, 
                minPi)) {
        double dualBoundValue, primalBoundValue, stepSize;
        vector <Edge> Eu(E), dualSolution, primalSolution;
        vector <double> G(S.size(), -1.0);
        (*totalIterations)++;
        /* Solving the Lagrangian problem with the current set of multipliers */
        for (unsigned int i = 0; i < S.size(); i++) {
            for (vector <Edge>::iterator it = Eu.begin(); it != Eu.end(); it++) {
                if ((it->u == S[i].e.u && it->v == S[i].e.u) || 
                        (it->u == S[i].f.u && it->v == S[i].f.v)) {
                    it->w += u[i];
                }
            }
        }
        dualBoundValue = kruskal(&dualSolution, n, Eu);
        for (vector <double>::iterator it = u.begin(); it != u.end(); it++) {
            dualBoundValue -= (*it);
        }
        if ((*bestDualBoundValue) < dualBoundValue) {
            (*bestDualBoundValue) = dualBoundValue;
            (*bestDualBoundIteration) = (*totalIterations);
            iterationsWithoutImprovment = 0;
        } else {
            iterationsWithoutImprovment++;
        }
        /* If bestDualBoundValue has not improved in the last N subgradient iterations */
        /* with the current value of π */
        if (iterationsWithoutImprovment >= N) {
            /* then halve π */
            pi /= 2.0;
            iterationsWithoutImprovment = 0;
        }
        /* Defining subgradients for the relaxed constraints, evaluated at the current solution */
        for (unsigned int i = 0; i < S.size(); i++) {
            for (vector <Edge>::iterator it = dualSolution.begin(); it != dualSolution.end(); 
                    it++) {
                if ((it->u == S[i].e.u && it->v == S[i].e.v) || 
                        (it->u == S[i].f.u && it->v == S[i].f.v)) {
                    G[i] += 1.0;
                }
            }
            if (fabs(u[i]) < EPSILON && G[i] < 0.0) {
                G[i] = 0.0;
            }
        }
        /* Obtaining a feasible primal solution from a (possibly unfeasible) dual solution */
        primalSolution = vector <Edge> (dualSolution);
        primalBoundValue = fixSolution(&primalSolution, n, E, S);
        if ((*bestPrimalBoundValue) > primalBoundValue) {
            (*bestPrimalBoundValue) = primalBoundValue;
            (*bestPrimalBoundIteration) = (*totalIterations);
            (*bestPrimalSolution) = primalSolution;
        }
        /* Defining a step size */
        stepSize = 0.0;
        for (vector <double>::iterator it = G.begin(); it != G.end(); it++) {
            stepSize += (*it) * (*it);
        }
        stepSize = (pi * ((*bestPrimalBoundValue) - dualBoundValue)) / stepSize;
        /* Updating Lagrange multipliers */
        for (unsigned int i = 0; i < u.size(); i++) {
            u[i] = max(0.0, u[i] + stepSize * G[i]);
        }
    }
    (*totalIterations)++;
    return true;
}

bool relaxLag2 (double * bestDualBoundValue, int * bestDualBoundIteration, int * totalIterations, 
        double * bestPrimalBoundValue, int * bestPrimalBoundIteration, 
        vector <Edge> * bestPrimalSolution, unsigned int n, vector <Edge> E, 
        vector <ConflictingPair> S, chrono :: high_resolution_clock :: time_point tBegin, 
        unsigned int timeLimit, double pi, unsigned int N, double minPi) {
    /* TODO */
    return true;
}

bool writeResult (char * resultFilePath, double bestDualBoundValue, int bestDualBoundIteration, 
        int totalIterations, double bestPrimalBoundValue, int bestPrimalBoundIteration, 
        vector <Edge> bestPrimalSolution) {
    ofstream resultFileStream(resultFilePath, ofstream::out);
    if (!resultFileStream.is_open()) {
        return false;
    }
    /* Uma linha com o melhor limitante dual encontrado impresso com 6 casas decimais. */
    resultFileStream << fixed << setprecision(6) << bestDualBoundValue << endl;
    /* Uma linha com a iteração do método do subgradiente (MS) */
    /* em que o melhor limitante dual */
    /* foi encontrado – considere que as iterações do MS são numeradas a partir de 0. */
    resultFileStream << bestDualBoundIteration << endl;
    /* Uma linha com a quantidade de iterações do MS executadas. */
    resultFileStream << totalIterations << endl;
    /* Uma linha com o melhor limitante primal encontrado. */
    resultFileStream << ((int) bestPrimalBoundValue) << endl;
    /* Uma linha com a iteração do MS em que o melhor */
    /* limitante primal foi encontrado – se o */
    /* melhor limitante primal encontrado foi calculado antes da execução do MS, */
    /* o número −1 deve ser impresso nesta linha. */
    resultFileStream << bestPrimalBoundIteration << endl;
    /* n−1 linhas, sendo n o número de vértices do grafo da instância de teste, */
    /* representando a árvore geradora correspondente ao melhor limitante primal encontrado, */
    /* cada linha identificando uma aresta e = (u, v) presente na árvore, */
    /* sendo 0 ≤ u, v ≤ n − 1, no seguinte formato: u v */
    for (vector <Edge>::iterator it = bestPrimalSolution.begin(); it != bestPrimalSolution.end(); 
            it++) {
        Edge e = *it;
        resultFileStream << e.u << ' ' << e.v << endl;
    }
    return true;
}

int main (int argc, char * argv[]) {
    chrono::high_resolution_clock::time_point tBegin = chrono::high_resolution_clock::now();
    if (argc != 4) {
        cerr << "Invalid arguments!" << endl;
        return 1;
    }
    unsigned int k = atoi(argv[1]);
    if (k < 1 || k > 2) {
        cerr << "Invalid relaxation number!" << endl;
        return 1;
    }

    unsigned int timeLimit, N;
    double pi, minPi;

    if (!readParameters(k, &timeLimit, &pi, &N, &minPi)) {
        cerr << "Error while reading parameters!" << endl;
        return 1;
    }

    unsigned int n;
    vector <Edge> E;
    vector <ConflictingPair> S;

    if (!readInput(&n, &E, &S, argv[2])) {
        cerr << "Error while reading input!" << endl;
    }

    double bestDualBoundValue, bestPrimalBoundValue;
    int bestDualBoundIteration, totalIterations, bestPrimalBoundIteration;
    vector <Edge> bestPrimalSolution;

    if (k == 1) {
        if (!relaxLag1(&bestDualBoundValue, &bestDualBoundIteration, &totalIterations, 
                    &bestPrimalBoundValue, &bestPrimalBoundIteration, &bestPrimalSolution, 
                    n, E, S, tBegin, timeLimit, pi, N, minPi)) {
            cerr << "Error while executing first Lagrangian Relaxation!" << endl;
            return 1;
        }
    } else {
        if (!relaxLag2(&bestDualBoundValue, &bestDualBoundIteration, &totalIterations, 
                    &bestPrimalBoundValue, &bestPrimalBoundIteration, &bestPrimalSolution, 
                    n, E, S, tBegin, timeLimit, pi, N, minPi)) {
            cerr << "Error while executing second Lagrangian Relaxation!" << endl;
            return 1;
        }
    }

    if (!writeResult(argv[3], bestDualBoundValue, bestDualBoundIteration, totalIterations, 
                bestPrimalBoundValue, bestPrimalBoundIteration, bestPrimalSolution)) {
        cerr << "Error while writing output!" << endl;
        return 1;
    }

    return 0;
}

