#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <utility>
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

bool isConnected (unsigned int n, vector <Edge> E) {
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

bool areEdgesExtremesEquals (Edge e, Edge f) {
    return (minmax(e.u, e.v) == minmax(f.u, f.v));
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
            if (areEdgesExtremesEquals((*S)[j].e, (*it))) {
                (*S)[j].e.w = it->w;
            }

            if (areEdgesExtremesEquals((*S)[j].f, (*it))) {
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
    if (fabs(bestDualBoundValue - bestPrimalBoundValue) < 1.0) {
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

    if (!isConnected(n, solution)) {
        return false;
    }

    for (vector <ConflictingPair>::iterator it = S.begin(); it != S.end(); it++) {
        unsigned int counter = 0;

        for (vector <Edge>::iterator it2 = solution.begin(); it2 != solution.end(); it2++) {
            if (areEdgesExtremesEquals(it->e, (*it2)) || areEdgesExtremesEquals(it->f, (*it2))) {
                counter++;
            }

            if (counter >= 2) {
                return false;
            }
        }
    }

    return true;
}

void removeEdge (vector <Edge> * E, Edge e) {
    vector <Edge>::iterator it = (*E).begin();
    while (it != (*E).end() && !areEdgesExtremesEquals((*it), e)) {
        it++;
    }
    if (it != (*E).end()) {
        (*E).erase(it);
    }
}

bool isBridge (unsigned int n, vector <Edge> E, Edge e) {
    Sets components;

    removeEdge (&E, e);
    components = connectedComponents(n, E);

    return !sameComponent(components, e.u, e.v);
}

unsigned int searchEdgeToRemove (unsigned int n, vector <Edge> E, vector <ConflictingPair> S, 
        vector <bool> isInPrimalSolution, vector <bool> isRemoved, vector <bool> isRemovable) {
    int result = E.size();
    vector <Edge> EnotRemoved;
    vector <unsigned int> conflictsCounter (E.size(), 0);
    unsigned int maxConflicts = 0, maxWeight = 0;

    for (unsigned int i = 0; i < E.size(); i++) {
        if (!isRemoved[i]) {
            EnotRemoved.push_back(E[i]);
        }
    }

    for (vector <ConflictingPair>::iterator it = S.begin(); it != S.end(); it++) {
        unsigned int counter = 0;
        for (unsigned int i = 0; i < E.size(); i++) {
            if (isInPrimalSolution[i]) {
                if (areEdgesExtremesEquals (it->e, E[i]) || areEdgesExtremesEquals(it->f, E[i])) {
                    counter++;
                }

                if (counter >= 2) {
                    conflictsCounter[i]++;
                    break;
                }
            }
        }
    }

    for (unsigned int i = 0; i < E.size(); i++) {
        if (isInPrimalSolution[i] && isRemovable[i] && conflictsCounter[i] > 0 && 
                (maxConflicts < conflictsCounter[i] || 
                 (maxConflicts == conflictsCounter[i] && maxWeight < E[i].w)) && 
                !isBridge(n, EnotRemoved, E[i])) {
            maxConflicts = conflictsCounter[i];
            maxWeight = E[i].w;
            result = i;
        }
    }

    return result;
}

unsigned int searchEdgeToInsert (unsigned int n, vector <Edge> E, vector <ConflictingPair> S, 
        vector <bool> isInPrimalSolution, vector <bool> isRemoved) {
    unsigned int result = E.size();
    vector <Edge> primalSolution;
    Sets components;

    for (unsigned int i = 0; i < E.size(); i++) {
        if (isInPrimalSolution[i]) {
            primalSolution.push_back(E[i]);
        }
    }

    components = connectedComponents(n, primalSolution);

    /* for each edge (u, v) ∈ E, taken in nondecreasing order by weight */
    for (unsigned int i = 0; i < E.size(); i++) {
        if (!isRemoved[i]) {
            /* Checking if the edge connects the solution */
            if (!sameComponent(components, E[i].u, E[i].v)) {
                /* Checking if the edge conflicts if another edge in solution */
                bool conflict = false;

                for (vector <ConflictingPair>::iterator it = S.begin(); 
                        it != S.end() && !conflict; it++) {
                    if (areEdgesExtremesEquals(it->e, E[i]) || 
                            areEdgesExtremesEquals(it->f, E[i])) {
                        for (vector <Edge>::iterator it2 = primalSolution.begin(); 
                                it2 != primalSolution.end(); it2++) {
                            if (areEdgesExtremesEquals(it->e, (*it2)) || 
                                    areEdgesExtremesEquals(it->f, (*it2))) {
                                conflict = true;
                            }
                        }
                    }
                }

                if (!conflict) {
                    result = i;
                    break;
                }
            }
        }
    }

    return result;
}

double fixSolution (vector <Edge> * primalSolution, unsigned int n, vector <Edge> E, 
        vector <ConflictingPair> S) {
    double primalBoundValue = 0.0;
    bool flag = false;
    vector <bool> isInPrimalSolution (E.size(), false), isRemoved (E.size(), false), 
           isRemovable (E.size(), true);

    /* Fixing the weight of solution's edges */
    for (unsigned int i = 0; i < E.size(); i++) {
        for (vector <Edge>::iterator it = (*primalSolution).begin(); 
                it != (*primalSolution).end(); it++) {
            if (areEdgesExtremesEquals(E[i], (*it))) {
                it->w = E[i].w;
                primalBoundValue += E[i].w;
            }
        }

        isRemovable[i] = !isBridge(n, E, E[i]);
    }

    while (!isFeasible(n, S, (*primalSolution)) && !flag) {
        unsigned int edgeToRemoveIndex;

        /* Converts the vector <Edge> codification into the vector <bool> codification */
        for (unsigned int i = 0; i < E.size(); i++) {
            for (vector <Edge>::iterator it = (*primalSolution).begin(); 
                    it != (*primalSolution).end(); it++) {
                if (areEdgesExtremesEquals(E[i], (*it))) {
                    isInPrimalSolution[i] = true;
                }
            }
        }

        /* Search for a conflicting edge to remove */
        edgeToRemoveIndex = searchEdgeToRemove(n, E, S, isInPrimalSolution, isRemoved, 
                isRemovable);

        /* If none found */
        if (edgeToRemoveIndex >= E.size()) {
            /* Set flag to stop loop */
            flag = true;
        } else {
            unsigned int edgeToInsertIndex;

            /* Remove edge from solution */
            isInPrimalSolution[edgeToRemoveIndex] = false;

            /* Remove edge from graph */
            isRemoved[edgeToRemoveIndex] = true;

            /* Search for a non-conflicting edge that re-connects the solution */
            edgeToInsertIndex = searchEdgeToInsert(n, E, S, isInPrimalSolution, isRemoved);

            /* If none found */
            if (edgeToInsertIndex >= E.size()) {
                /* Re-add edge to solution */
                isInPrimalSolution[edgeToRemoveIndex] = true;

                /* Re-add edge to graph */
                isRemoved[edgeToRemoveIndex] = false;

                /* Sets edge as non removable */
                isRemovable[edgeToRemoveIndex] = false;
            } else {
                /* Add edge to solution */
                isInPrimalSolution[edgeToInsertIndex] = true;
            }
        }

        /* Converts the vector <bool> codification into the vector <Edge> codification */
        (*primalSolution) = vector <Edge> ();
        primalBoundValue = 0.0;
        for (unsigned int i = 0; i < E.size(); i++) {
            if (isInPrimalSolution[i]) {
                (*primalSolution).push_back(E[i]);
                primalBoundValue += E[i].w;
            }
        }
    }

    return primalBoundValue;
}

void writeOutput (int iteration, chrono::high_resolution_clock::time_point tBegin, 
        double dualBoundValue, double primalBoundValue) {
    chrono::high_resolution_clock::time_point tCurrent = chrono::high_resolution_clock::now();
    chrono::seconds elapsedTime = chrono::duration_cast <chrono::seconds> (tCurrent - tBegin);

    cout << iteration << " ";
    cout << elapsedTime.count() << " ";
    cout << fixed << setprecision(6) << dualBoundValue << " ";

    if (!isnan(primalBoundValue)) {
        cout << ((int) primalBoundValue) << endl;
    } else {
        cout << "NaN" << endl;
    }
}

bool relaxLag1 (double * bestDualBoundValue, int * bestDualBoundIteration, int * totalIterations, 
        double * bestPrimalBoundValue, int * bestPrimalBoundIteration, 
        vector <Edge> * bestPrimalSolution, unsigned int n, vector <Edge> E, 
        vector <ConflictingPair> S, chrono::high_resolution_clock::time_point tBegin, 
        unsigned int timeLimit, double pi, unsigned int N, double minPi) {
    vector <double> u;
    unsigned int iterationsWithoutImprovment;
    double dualBoundValue, primalBoundValue;
    vector <Edge> dualSolution, primalSolution;

    (*bestDualBoundValue) = -INFINITE;
    (*bestDualBoundIteration) = -1;
    (*totalIterations) = 0;
    (*bestPrimalBoundValue) = INFINITE;
    (*bestPrimalBoundIteration) = -1;

    /* sort the edges of E into nondecreasing order by weight w */
    sort(E.begin(), E.end(), comparator);

    dualBoundValue = kruskal(&dualSolution, n, E);

    if ((*bestDualBoundValue) < dualBoundValue) {
        (*bestDualBoundValue) = dualBoundValue;
    }

    primalSolution = vector <Edge> (dualSolution);
    primalBoundValue = fixSolution(&primalSolution, n, E, S);

    if (isFeasible(n, S, primalSolution) && (*bestPrimalBoundValue) > primalBoundValue) {
        (*bestPrimalBoundValue) = primalBoundValue;
        (*bestPrimalSolution) = primalSolution;
    }

    if (isFeasible(n, S, primalSolution)) {
        writeOutput((*totalIterations), tBegin, dualBoundValue, primalBoundValue);
    } else {
        writeOutput((*totalIterations), tBegin, dualBoundValue, nan(""));
    }

    u = initialLagrangeMultipliers(S.size());
    iterationsWithoutImprovment = 0;

    while (!termination(tBegin, timeLimit, (*bestDualBoundValue), (*bestPrimalBoundValue), pi, 
                minPi)) {
        double stepSize;
        vector <Edge> Eu(E);
        vector <double> G(S.size(), -1.0);

        (*totalIterations)++;

        /* Solving the Lagrangian problem with the current set of multipliers */
        for (unsigned int i = 0; i < S.size(); i++) {
            for (vector <Edge>::iterator it = Eu.begin(); it != Eu.end(); it++) {
                if (areEdgesExtremesEquals((*it), S[i].e) || 
                        areEdgesExtremesEquals((*it), S[i].f)) {
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
                if (areEdgesExtremesEquals((*it), S[i].e) || 
                        areEdgesExtremesEquals((*it), S[i].f)) {
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

        if (isFeasible(n, S, primalSolution) && (*bestPrimalBoundValue) > primalBoundValue) {
            (*bestPrimalBoundValue) = primalBoundValue;
            (*bestPrimalBoundIteration) = (*totalIterations);
            (*bestPrimalSolution) = primalSolution;
        }

        if (isFeasible(n, S, primalSolution)) {
            writeOutput((*totalIterations), tBegin, dualBoundValue, primalBoundValue);
        } else {
            writeOutput((*totalIterations), tBegin, dualBoundValue, nan(""));
        }

        /* Defining a step size */
        stepSize = 0.0;

        for (vector <double>::iterator it = G.begin(); it != G.end(); it++) {
            stepSize += (*it) * (*it);
        }

        stepSize = (pi * ((1.05 * (*bestPrimalBoundValue)) - dualBoundValue)) / stepSize;

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
    vector < vector <ConflictingPair> > Se;
    unsigned int eStarIndex;
    vector <ConflictingPair> SminusSeStar;
    vector <Edge> EeStar, EminusEeStar, dualSolution, primalSolution;
    vector <double> u;
    unsigned int iterationsWithoutImprovment;
    double dualBoundValue, primalBoundValue;

    (*bestDualBoundValue) = -INFINITE;
    (*bestDualBoundIteration) = -1;
    (*totalIterations) = 0;
    (*bestPrimalBoundValue) = INFINITE;
    (*bestPrimalBoundIteration) = -1;

    /* sort the edges of E into nondecreasing order by weight w */
    sort(E.begin(), E.end(), comparator);

    dualBoundValue = kruskal(&dualSolution, n, E);

    if ((*bestDualBoundValue) < dualBoundValue) {
        (*bestDualBoundValue) = dualBoundValue;
    }

    primalSolution = vector <Edge> (dualSolution);
    primalBoundValue = fixSolution(&primalSolution, n, E, S);

    if (isFeasible(n, S, primalSolution) && (*bestPrimalBoundValue) > primalBoundValue) {
        (*bestPrimalBoundValue) = primalBoundValue;
        (*bestPrimalSolution) = primalSolution;
    }

    if (isFeasible(n, S, primalSolution)) {
        writeOutput((*totalIterations), tBegin, dualBoundValue, primalBoundValue);
    } else {
        writeOutput((*totalIterations), tBegin, dualBoundValue, nan(""));
    }

    /* Denote por Se o conjunto de pares conflitantes de S envolvendo a aresta e */
    Se = vector < vector <ConflictingPair> > (E.size());

    for (unsigned int i = 0; i < E.size(); i++) {
        for (vector <ConflictingPair>::iterator it = S.begin(); it != S.end(); it++) {
            if (areEdgesExtremesEquals(E[i], it->e) || areEdgesExtremesEquals(E[i], it->f)) {
                Se[i].push_back((*it));
            }
        }
    }

    /* Seja ainda e* a aresta de E com menor custo para a qual Se não é vazio */
    for (eStarIndex = 0; eStarIndex < E.size(); eStarIndex++) {
        if (Se[eStarIndex].size() > 0) {
            break;
        }
    }

    for (vector <ConflictingPair>::iterator it = S.begin(); it != S.end(); it++) {
        bool found = false;

        for (vector <ConflictingPair>::iterator it2 = Se[eStarIndex].begin(); 
                it2 != Se[eStarIndex].end() && !found; it2++) {
            if (areEdgesExtremesEquals(it->e, it2->e) && areEdgesExtremesEquals(it->f, it2->f)) {
                found = true;
            }
        }

        if (!found) {
            SminusSeStar.push_back((*it));
        }
    }

    for (vector <Edge>::iterator it = E.begin(); it != E.end(); it++) {
        if (!areEdgesExtremesEquals((*it), E[eStarIndex])) {
            bool conflictsWithEStar = false;

            for (vector <ConflictingPair>::iterator it2 = Se[eStarIndex].begin(); 
                    it2 != Se[eStarIndex].end() && !conflictsWithEStar; it2++) {
                if (areEdgesExtremesEquals((*it), it2->e) || 
                        areEdgesExtremesEquals((*it), it2->f)) {
                    conflictsWithEStar = true;
                }
            }

            if (conflictsWithEStar) {
                EeStar.push_back((*it));
            }
        }
    }

    for (vector <Edge>::iterator it = E.begin(); it != E.end(); it++) {
        bool found = false;

        for (vector <Edge>::iterator it2 = EeStar.begin(); it2 != EeStar.end() && !found; it2++) {
            if (areEdgesExtremesEquals((*it), (*it2))) {
                found = true;
            }
        }

        if (!found) {
            EminusEeStar.push_back((*it));
        }
    }

    u = initialLagrangeMultipliers(SminusSeStar.size());
    iterationsWithoutImprovment = 0;

    while (!termination(tBegin, timeLimit, (*bestDualBoundValue), (*bestPrimalBoundValue), pi, 
                minPi)) {
        double stepSize;
        vector <Edge> Eu(EminusEeStar);
        vector <double> G(SminusSeStar.size(), -1.0);

        (*totalIterations)++;

        /* Solving the Lagrangian problem with the current set of multipliers */
        for (unsigned int i = 0; i < SminusSeStar.size(); i++) {
            for (vector <Edge>::iterator it = Eu.begin(); it != Eu.end(); it++) {
                if (areEdgesExtremesEquals((*it), SminusSeStar[i].e) || 
                        areEdgesExtremesEquals((*it), SminusSeStar[i].f)) {
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
        for (unsigned int i = 0; i < SminusSeStar.size(); i++) {
            for (vector <Edge>::iterator it = dualSolution.begin(); it != dualSolution.end(); 
                    it++) {
                if (areEdgesExtremesEquals((*it), SminusSeStar[i].e) || 
                        areEdgesExtremesEquals((*it), SminusSeStar[i].f)) {
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

        if (isFeasible(n, S, primalSolution) && (*bestPrimalBoundValue) > primalBoundValue) {
            (*bestPrimalBoundValue) = primalBoundValue;
            (*bestPrimalBoundIteration) = (*totalIterations);
            (*bestPrimalSolution) = primalSolution;
        }

        if (isFeasible(n, S, primalSolution)) {
            writeOutput((*totalIterations), tBegin, dualBoundValue, primalBoundValue);
        } else {
            writeOutput((*totalIterations), tBegin, dualBoundValue, nan(""));
        }

        /* Defining a step size */
        stepSize = 0.0;

        for (vector <double>::iterator it = G.begin(); it != G.end(); it++) {
            stepSize += (*it) * (*it);
        }
        
        stepSize = (pi * ((1.05 * (*bestPrimalBoundValue)) - dualBoundValue)) / stepSize;

        /* Updating Lagrange multipliers */
        for (unsigned int i = 0; i < u.size(); i++) {
            u[i] = max(0.0, u[i] + stepSize * G[i]);
        }
    }

    (*totalIterations)++;

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

vector < list < pair <unsigned int, double> > > toAdjacencyList (unsigned int n, 
        vector <Edge> E) {
    vector < list < pair <unsigned int, double> > > adj (n);

    for (vector <Edge>::iterator it = E.begin(); it != E.end(); it++) {
        Edge e = (*it);
        adj[e.u].push_back(make_pair(e.v, e.w));
        adj[e.v].push_back(make_pair(e.u, e.w));
    }

    return adj;
}

void removeEdgeFromAdjacencyList(vector < list < pair <unsigned int, double> > > adj, Edge e) {
    /* Searchs for edge (u, v) */
    list < pair <unsigned int, double> >::iterator it = adj[e.u].begin();

    while (it != adj[e.u].end() && it->first != e.v) {
        it++;
    }

    /* If found, remove edge (u, v) */
    if (it != adj[e.u].end()) {
        adj[e.u].erase(it);
    }

    /* Searched for edge (v, u) */
    it = adj[e.v].begin();

    while (it != adj[e.v].end() && it->first != e.u) {
        it++;
    }

    /* If found, remove edge (v, u) */
    if (it != adj[e.v].end()) {
        adj[e.v].erase(it);
    }
}

void bridgesAux (vector <list < pair <unsigned int, double> > > adj, unsigned int * time, 
        vector <bool> * visited, vector <unsigned int> * disc, vector <unsigned int> * low, 
        vector <unsigned int> * parent, unsigned int u, vector <Edge> * bridges) {
    (*visited)[u] = true;
    (*time)++;
    (*disc)[u] = (*time);
    (*low)[u] = (*time);

    for (list < pair <unsigned int, double> >::iterator it = adj[u].begin(); it != adj[u].end(); 
            it++) {
        unsigned int v = (*it).first;
        double w = (*it).second;

        if (!(*visited)[v]) {
            (*parent)[v] = u;
            bridgesAux(adj, time, visited, disc, low, parent, v, bridges);

            if ((*low)[u] > (*low)[v]) {
                (*low)[u] = (*low)[v];
            }

            if ((*low)[v] > (*disc)[u]) {
                Edge e;
                e.u = u;
                e.v = v;
                e.w = w;
                (*bridges).push_back(e);
            }
        } else if (v != (*parent)[u] && (*low)[u] > (*disc)[v]) {
            (*low)[u] = (*disc)[v];
        }
    }
}

vector <Edge> getBridges (vector <list < pair <unsigned int, double> > > adj) {
    vector <Edge> result;
    unsigned int time = 0;
    vector <bool> visited (adj.size(), false);
    vector <unsigned int> disc (adj.size(), 0), low (adj.size(), 0);
    vector <unsigned int> parent (adj.size(), adj.size());

    for (unsigned int u = 0; u < adj.size(); u++) {
        if (!visited[u]) {
            bridgesAux(adj, &time, &visited, &disc, &low, &parent, u, &result);
        }
    }

    return result;
}

bool isConflictingWithSomeEdge (vector <ConflictingPair> S, Edge e) {
    for (vector <ConflictingPair>::iterator it = S.begin(); it != S.end(); it++) {
        if (areEdgesExtremesEquals(e, it->e) || areEdgesExtremesEquals(e, it->f)) {
            return true;
        }
    }

    return false;
}

bool areEdgesConflicting (vector <ConflictingPair> S, Edge e, Edge f) {
    for (vector <ConflictingPair>::iterator it = S.begin(); it != S.end(); it++) {
        if ((areEdgesExtremesEquals(e, it->e) && areEdgesExtremesEquals(f, it->f)) || 
                    (areEdgesExtremesEquals(e, it->f) && areEdgesExtremesEquals(f, it->e))) {
            return true;
        }
    }

    return false;
}

void preProcessingPhase1 (unsigned int n, vector <Edge> * E, vector <ConflictingPair> * S, 
        vector <Edge> * fixedEdges) {
    bool flag = false;

    vector < list < pair <unsigned int, double> > > adj = toAdjacencyList (n, (*E));

    while (!flag && isConnected(n, (*E))) {
        vector <Edge> bridges = getBridges(adj);

        for (vector <Edge>::iterator it = (*fixedEdges).begin(); it != (*fixedEdges).end(); it++) {
            removeEdge(&bridges, (*it));
        }

        /* For each bridge e not yet fixed */
        for (vector <Edge>::iterator it = bridges.begin(); it != bridges.end(); it++) {
            Edge e = (*it);

            /* Add e to the fixedEdges */
            (*fixedEdges).push_back(e);

            /* For each conflicting pair containing e */
            for (vector <ConflictingPair>::iterator it2 = (*S).begin(); it2 != (*S).end();) {
                if (areEdgesExtremesEquals(e, it2->e)) {
                    /* Remove the conflicing edge from the graph */
                    removeEdge(E, it2->f);
                    removeEdgeFromAdjacencyList(adj, it2->f);

                    /* Remove the conflicting pair */
                    it2 = (*S).erase(it2);
                } else if (areEdgesExtremesEquals(e, it2->f)) {
                    /* Remove the conflicing edge from the graph */
                    removeEdge(E, it2->e);
                    removeEdgeFromAdjacencyList(adj, it2->e);

                    /* Remove the conflicting pair */
                    it2 = (*S).erase(it2);
                } else {
                    it2++;
                }
            }
        }

        if (bridges.size() <= 0) {
            flag = true;
        }
    }
}

bool preProcessingPhase2 (unsigned int n, vector <Edge> * E, vector <ConflictingPair> * S, 
        vector <Edge> * fixedEdges) {
    bool result = false;

    /* For each edge e that is conflicting with some edge */
    for (vector <Edge>::iterator it = (*E).begin(); it != (*E).end();) {
        Edge e = (*it);

        if (isConflictingWithSomeEdge((*S), e)) {
            /* Create a copy of the graph */
            vector <Edge> Eprime ((*E));

            /* Create a copy of the conflicting pairs */
            vector <ConflictingPair> Sprime ((*S));

            /* Create a copy of the fixed edges */
            vector <Edge> fixedEdgesPrime ((*fixedEdges));

            /* Fix e in the copied instance */
            /* by removing the edges that conflict with e from the copied graph */
            /* and removing conflicts involving e from the copied conflicting pairs */
            for (vector <ConflictingPair>::iterator it2 = Sprime.begin(); it2 != Sprime.end();) {
                if (areEdgesExtremesEquals(e, it2->e)) {
                    removeEdge(&Eprime, it2->f);

                    it2 = Sprime.erase(it2);
                } else if (areEdgesExtremesEquals(e, it2->f)) {
                    removeEdge(&Eprime, it2->e);

                    it2 = Sprime.erase(it2);
                } else {
                    it2++;
                }
            }

            /* Apply pre-processing phase 1 in the copied instance */
            /* in order to fix the new bridges */
            preProcessingPhase1(n, &Eprime, &Sprime, &fixedEdgesPrime);

            /* If the copied graph becomes disconnected */
            if (!isConnected(n, Eprime)) {
                /* Removes conflicts involving e */
                for (vector <ConflictingPair>::iterator it2 = (*S).begin(); it2 != (*S).end();) {
                    if (areEdgesExtremesEquals(e, it2->e) || areEdgesExtremesEquals(e, it2->f)) {
                        it2 = (*S).erase(it2);
                    } else {
                        it2++;
                    }
                }

                /* Remove e from the graph */
                it = (*E).erase(it);

                /* Update return value */
                result = true;
            } else {
                it++;
            }
        } else {
            it++;
        }
    }

    return result;
}

bool preProcessingPhase3 (unsigned int n, vector <Edge> * E, vector <ConflictingPair> * S, 
        vector <Edge> * fixedEdges) {
    bool result = false;
    /* TODO */
    return result;
}

void preProcessing (unsigned int n, vector <Edge> * E, vector <ConflictingPair> * S, 
        vector <Edge> * fixedEdges) {
    do {
        do {
            preProcessingPhase1(n, E, S, fixedEdges);
        } while (isConnected(n, (*E)) && preProcessingPhase2(n, E, S, fixedEdges));
    } while (isConnected(n, (*E)) && preProcessingPhase3(n, E, S, fixedEdges));
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

    vector <Edge> fixedEdges;
    preProcessing (n, &E, &S, &fixedEdges);

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

