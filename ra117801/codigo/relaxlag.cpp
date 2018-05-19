#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
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
    sort(E.begin(), E.end(), edgeWeightComparator);

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

bool areEdgesExtremesEquals (Edge e, Edge f) {
    return (minmax(e.u, e.v) == minmax(f.u, f.v));
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

bool readParameters (unsigned int k, unsigned int * totalTimeLimit, 
        unsigned int * preProcessingTimeLimit, unsigned int * fixSolutionTimeLimit, double * pi, 
        unsigned int * N, double * minPi) {
    ifstream parameterFileStream(k == 1 ? "param1" : "param2", ifstream::in);

    if (!parameterFileStream.is_open()) {
        return false;
    }

    /* Total time limit in seconds */
    parameterFileStream >> (*totalTimeLimit);

    /* Pre-processing time limit in seconds */
    parameterFileStream >> (*preProcessingTimeLimit);

    /* Fix solution time limit in seconds */
    parameterFileStream >> (*fixSolutionTimeLimit);

    /* π */
    parameterFileStream >> (*pi);

    /* N */
    parameterFileStream >> (*N);

    /* minPi */
    parameterFileStream >> (*minPi);

    return true;
}

bool readInput (unsigned int * n, vector <Edge> * E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> * S, 
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

    bool (*conflictingPairComparatorPointer) (ConflictingPair, ConflictingPair) = 
        conflictingPairComparator;
    (*S) = set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> (conflictingPairComparatorPointer);

    /* m linhas, cada uma identificando uma aresta e = (u, v), de custo ce, sendo 0 ≤ u, */
    /* v ≤ n − 1 e 0 ≤ ce ≤ 500, no seguinte formato: u v ce */
    for (unsigned int j = 0; j < m; j++) {
        inputFileStream >> (*E)[j].u >> (*E)[j].v >> (*E)[j].w;
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

        for (vector <Edge>::iterator it = (*E).begin(); it != (*E).end(); it++) {
            if (areEdgesExtremesEquals(e, (*it))) {
                e.w = it->w;
            }

            if (areEdgesExtremesEquals(f, (*it))) {
                f.w = it->w;
            }
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

vector <double> initialLagrangeMultipliers (int m) {
    vector <double> result (m, 1.0);
    return result;
}

bool timeLimitExceeded (chrono::high_resolution_clock::time_point tBegin, 
        unsigned int timeLimit) {
    chrono::high_resolution_clock::time_point tCurrent = chrono::high_resolution_clock::now();
    chrono::seconds elapsedTime = chrono::duration_cast <chrono::seconds> (tCurrent - tBegin);

    if ((unsigned int) elapsedTime.count() >= timeLimit) {
        return true;
    }

    return false;
}

bool termination (chrono::high_resolution_clock::time_point tBegin, unsigned int totalTimeLimit, 
        double bestDualBoundValue, double bestPrimalBoundValue, double pi, double minPi) {
    if (timeLimitExceeded(tBegin, totalTimeLimit)) {
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

bool isFeasible (unsigned int n, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        vector <Edge> solution) {
    if (solution.size() != n - 1) {
        return false;
    }

    if (!isConnected(n, solution)) {
        return false;
    }

    for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it = 
            S.begin(); it != S.end(); it++) {
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

bool searchEdgeToRemove (Edge * edgeToRemove, unsigned int n, vector <Edge> E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        set <Edge, bool (*) (Edge, Edge)> fixedEdges, vector <Edge> primalSolution) {
    bool result = false;

    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    set <Edge, bool (*) (Edge, Edge)> removableEdges (edgeComparatorPointer);

    bool (*conflictingPairComparatorPointer) (ConflictingPair, ConflictingPair) = 
        conflictingPairComparator;
    set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> Sative 
        (conflictingPairComparatorPointer);

    for (vector <Edge>::iterator it = primalSolution.begin(); it != primalSolution.end(); it++) {
        Edge e = (*it);
        for (vector <Edge>::iterator it2 = it + 1; it2 != primalSolution.end(); it2++) {
            Edge f = (*it2);

            ConflictingPair cp;

            if (edgeComparator(e, f)) {
                cp.e = e;
                cp.f = f;
            } else {
                cp.e = f;
                cp.f = e;
            }

            if (S.find(cp) != S.end()) {
                if (fixedEdges.find(e) == fixedEdges.end()) {
                    removableEdges.insert(e);
                }

                if (fixedEdges.find(f) == fixedEdges.end()) {
                    removableEdges.insert(f);
                }

                Sative.insert(cp);
            }
        }
    }

    unsigned int maxConflictsCounter = 0;
    double maxWeight = -INFINITE;

    for (set <Edge, bool (*) (Edge, Edge)>::iterator it = removableEdges.begin(); 
            it != removableEdges.end(); it++) {
        Edge e = (*it);
        unsigned int conflictsCounter = 0;

        for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it2 = 
                Sative.begin(); it2 != Sative.end(); it2++) {
            if (areEdgesExtremesEquals(e, it2->e) || areEdgesExtremesEquals(e, it2->f)) {
                conflictsCounter++;
            }
        }

        if (maxConflictsCounter < conflictsCounter || 
                (maxConflictsCounter == conflictsCounter && maxWeight < e.w)) {
            maxConflictsCounter = conflictsCounter;
            maxWeight = e.w;
            (*edgeToRemove) = e;
            result = true;
        }
    }

    return result;
}

bool searchEdgeToInsert (Edge * edgeToInsert, unsigned int n, vector <Edge> E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        vector <Edge> primalSolution, set <Edge, bool (*) (Edge, Edge)> removedEdges) {
    bool result = false;
    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    set <Edge, bool (*) (Edge, Edge)> insertableEdges (E.begin(), E.end(), edgeComparatorPointer);

    for (set <Edge, bool (*) (Edge, Edge)>::iterator it = removedEdges.begin(); 
            it != removedEdges.end(); it++) {
        Edge e = (*it);

        insertableEdges.erase(e);
    }

    for (vector <Edge>::iterator it = primalSolution.begin(); it != primalSolution.end(); it++) {
        Edge e = (*it);

        insertableEdges.erase(e);

        for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it2 = 
                S.begin(); it2 != S.end(); it2++) {
            if (areEdgesExtremesEquals(e, it2->e)) {
                insertableEdges.erase(it2->f);
            } else if (areEdgesExtremesEquals(e, it2->f)) {
                insertableEdges.erase(it2->e);
            }
        }
    }

    if (insertableEdges.size() > 0) {
        Sets components = connectedComponents(n, primalSolution);
        double minWeight = INFINITE;

        for (set <Edge, bool (*) (Edge, Edge)>::iterator it = insertableEdges.begin(); 
                it != insertableEdges.end(); it++) {
            Edge e = (*it);

            if (!sameComponent(components, e.u, e.v) && minWeight > e.w) {
                minWeight = e.w;
                (*edgeToInsert) = e;
                result = true;
            }
        }
    }

    return result;
}

double fixSolution (vector <Edge> * primalSolution, unsigned int n, vector <Edge> E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        set <Edge, bool (*) (Edge, Edge)> fixedEdges, unsigned int fixSolutionTimeLimit) {
    chrono::high_resolution_clock::time_point tBeginFixSolution =
        chrono::high_resolution_clock::now();
    double primalBoundValue = 0.0;

    for (vector <Edge>::iterator it = (*primalSolution).begin(); it != (*primalSolution).end(); 
            it++) {
        for (vector <Edge>::iterator it2 = E.begin(); it2 != E.end(); it2++) {
            if (areEdgesExtremesEquals((*it), (*it2))) {
                it->w = it2->w;
                primalBoundValue += it->w;
            }
        }
    }

    Edge edgeToRemove;
    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    set <Edge, bool (*) (Edge, Edge)> removedEdges (edgeComparatorPointer);

    while (!timeLimitExceeded(tBeginFixSolution, fixSolutionTimeLimit) && 
            searchEdgeToRemove (&edgeToRemove, n, E, S, fixedEdges, (*primalSolution))) {
        Edge edgeToInsert;

        primalBoundValue -= edgeToRemove.w;
        removeEdge(primalSolution, edgeToRemove);
        removedEdges.insert(edgeToRemove);

        if (searchEdgeToInsert (&edgeToInsert, n, E, S, (*primalSolution), removedEdges)) {
            primalBoundValue += edgeToInsert.w;
            (*primalSolution).push_back(edgeToInsert);
        } else {
            primalBoundValue += edgeToRemove.w;
            (*primalSolution).push_back(edgeToRemove);

            break;
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

double constructiveHeuristic (vector <Edge> * primalSolution, unsigned int n, vector <Edge> E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S) {
    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    set <Edge, bool (*) (Edge, Edge)> conflictingEdges (edgeComparatorPointer);

    for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it = 
            S.begin(); it != S.end(); it++) {
        conflictingEdges.insert(it->e);
        conflictingEdges.insert(it->f);
    }

    /* While S = ∅ */
    while (S.size() > 0) {
        unsigned int maxConflictsCounter = 0;
        Edge maxConflictingEdge;

        /* In the graph G we choose an e ∈ E */
        /* which appears in the maximum number of conflict pairs */
        for (set <Edge>::iterator it = conflictingEdges.begin(); it != conflictingEdges.end();
                it++) {
            Edge e = (*it);
            unsigned int conflictsCounter = 0;

            for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it2 
                    = S.begin(); it2 != S.end(); it2++) {
                if (areEdgesExtremesEquals(e, it2->e) || areEdgesExtremesEquals(e, it2->f)) {
                    conflictsCounter++;
                }
            }

            if (maxConflictsCounter < conflictsCounter) {
                maxConflictsCounter = conflictsCounter;
                maxConflictingEdge = (*it);
            }
        }

        /* Delete e from G */
        removeEdge(&E, maxConflictingEdge);

        /* And delete from S all the conflict pairs containing e */
        for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it = 
                S.begin(); it != S.end();) {
            if (areEdgesExtremesEquals(maxConflictingEdge, it->e) ||
                    areEdgesExtremesEquals(maxConflictingEdge, it->f)) {
                it = S.erase(it);
            } else {
                it++;
            }
        }
    }

    /* If G is still connected */
    if (isConnected(n, E)) {
        /* then the cost of its minimum spanning tree is an upper bound to the MSTC */
        return kruskal (primalSolution, n, E);
    } else {
        return 0;
    }
}

bool relaxLag1 (double * bestDualBoundValue, int * bestDualBoundIteration, int * totalIterations, 
        double * bestPrimalBoundValue, int * bestPrimalBoundIteration, 
        vector <Edge> * bestPrimalSolution, unsigned int n, vector <Edge> E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        set <Edge, bool (*) (Edge, Edge)> fixedEdges, 
        chrono::high_resolution_clock::time_point tBegin, unsigned int totalTimeLimit, 
        unsigned int fixSolutionTimeLimit, double pi, unsigned int N, double minPi) {
    vector <double> u;
    unsigned int iterationsWithoutImprovment;
    double dualBoundValue, primalBoundValue;
    vector <Edge> dualSolution, primalSolution;

    /* sort the edges of E into nondecreasing order by weight w */
    sort(E.begin(), E.end(), edgeWeightComparator);

    dualBoundValue = kruskal(&dualSolution, n, E);

    if ((*bestDualBoundValue) < dualBoundValue) {
        (*bestDualBoundValue) = dualBoundValue;
    }

    primalSolution = vector <Edge> (dualSolution);

    primalBoundValue = fixSolution(&primalSolution, n, E, S, fixedEdges, fixSolutionTimeLimit);

    if (isFeasible(n, S, primalSolution) && (*bestPrimalBoundValue) > primalBoundValue) {
        (*bestPrimalBoundValue) = primalBoundValue;
        (*bestPrimalSolution) = primalSolution;
    }

    if (isFeasible(n, S, primalSolution)) {
        writeOutput((*totalIterations), tBegin, dualBoundValue, primalBoundValue);
    } else {
        writeOutput((*totalIterations), tBegin, dualBoundValue, nan(""));
    }

    constructiveHeuristic(&primalSolution, n, E, S);

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

    while (!termination(tBegin, totalTimeLimit, (*bestDualBoundValue), (*bestPrimalBoundValue), 
                pi, minPi)) {
        double stepSize;
        vector <Edge> Eu(E);
        vector <double> G(S.size(), -1.0);

        /* Solving the Lagrangian problem with the current set of multipliers */
        unsigned int i = 0;
        for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it = 
                S.begin(); it != S.end(); it++, i++) {
            for (vector <Edge>::iterator it2 = Eu.begin(); it2 != Eu.end(); it2++) {
                if (areEdgesExtremesEquals((*it2), it->e) || 
                        areEdgesExtremesEquals((*it2), it->f)) {
                    it2->w += u[i];
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
        i = 0;
        for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it = 
                S.begin(); it != S.end(); it++, i++) {
            for (vector <Edge>::iterator it2 = dualSolution.begin(); it2 != dualSolution.end(); 
                    it2++) {
                if (areEdgesExtremesEquals((*it2), it->e) || 
                        areEdgesExtremesEquals((*it2), it->f)) {
                    G[i] += 1.0;
                }
            }

            if (fabs(u[i]) < EPSILON && G[i] < 0.0) {
                G[i] = 0.0;
            }
        }

        /* Obtaining a feasible primal solution from a (possibly unfeasible) dual solution */
        primalSolution = vector <Edge> (dualSolution);
        primalBoundValue = fixSolution(&primalSolution, n, E, S, fixedEdges, 
                fixSolutionTimeLimit);

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

        (*totalIterations)++;
    }

    (*totalIterations)++;

    return true;
}

bool relaxLag2 (double * bestDualBoundValue, int * bestDualBoundIteration, int * totalIterations, 
        double * bestPrimalBoundValue, int * bestPrimalBoundIteration, 
        vector <Edge> * bestPrimalSolution, unsigned int n, vector <Edge> E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        set <Edge, bool (*) (Edge, Edge)> fixedEdges, 
        chrono::high_resolution_clock::time_point tBegin, unsigned int totalTimeLimit, 
        unsigned int fixSolutionTimeLimit, double pi, unsigned int N, double minPi) {
    vector < set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> > Se;
    unsigned int eStarIndex;
    set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> SminusSeStar (S);
    vector <Edge> EeStar, EminusEeStar, dualSolution, primalSolution;
    vector <double> u;
    unsigned int iterationsWithoutImprovment;
    double dualBoundValue, primalBoundValue;

    /* sort the edges of E into nondecreasing order by weight w */
    sort(E.begin(), E.end(), edgeWeightComparator);

    dualBoundValue = kruskal(&dualSolution, n, E);

    if ((*bestDualBoundValue) < dualBoundValue) {
        (*bestDualBoundValue) = dualBoundValue;
    }

    primalSolution = vector <Edge> (dualSolution);
    primalBoundValue = fixSolution(&primalSolution, n, E, S, fixedEdges, fixSolutionTimeLimit);

    if (isFeasible(n, S, primalSolution) && (*bestPrimalBoundValue) > primalBoundValue) {
        (*bestPrimalBoundValue) = primalBoundValue;
        (*bestPrimalSolution) = primalSolution;
    }

    if (isFeasible(n, S, primalSolution)) {
        writeOutput((*totalIterations), tBegin, dualBoundValue, primalBoundValue);
    } else {
        writeOutput((*totalIterations), tBegin, dualBoundValue, nan(""));
    }

    constructiveHeuristic(&primalSolution, n, E, S);

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
    Se = vector < set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> > (E.size());

    for (unsigned int i = 0; i < E.size(); i++) {
        for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it = 
                S.begin(); it != S.end(); it++) {
            if (areEdgesExtremesEquals(E[i], it->e) || areEdgesExtremesEquals(E[i], it->f)) {
                Se[i].insert((*it));
            }
        }
    }

    /* Seja ainda e* a aresta de E com menor custo para a qual Se não é vazio */
    for (eStarIndex = 0; eStarIndex < E.size(); eStarIndex++) {
        if (Se[eStarIndex].size() > 0) {
            break;
        }
    }

    for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it = 
            Se[eStarIndex].begin(); it != Se[eStarIndex].end(); it++) {
        SminusSeStar.erase((*it));
    }

    for (vector <Edge>::iterator it = E.begin(); it != E.end(); it++) {
        Edge e = (*it);

        if (!areEdgesExtremesEquals(e, E[eStarIndex])) {
            ConflictingPair cp;

            if (edgeComparator(e, E[eStarIndex])) {
                cp.e = e;
                cp.f = E[eStarIndex];
            } else {
                cp.e = E[eStarIndex];
                cp.f = e;
            }

            if (Se[eStarIndex].find(cp) != Se[eStarIndex].end()) {
                EeStar.push_back(e);
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

    while (!termination(tBegin, totalTimeLimit, (*bestDualBoundValue), (*bestPrimalBoundValue), 
                pi, minPi)) {
        double stepSize;
        vector <Edge> Eu(EminusEeStar);
        vector <double> G(SminusSeStar.size(), -1.0);

        /* Solving the Lagrangian problem with the current set of multipliers */
        unsigned int i = 0;
        for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it = 
                SminusSeStar.begin(); it != SminusSeStar.end(); it++, i++) {
            for (vector <Edge>::iterator it2 = Eu.begin(); it2 != Eu.end(); it2++) {
                if (areEdgesExtremesEquals((*it2), it->e) || 
                        areEdgesExtremesEquals((*it2), it->f)) {
                    it2->w += u[i];
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
        i = 0;
        for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it = 
                SminusSeStar.begin(); it != SminusSeStar.end(); it++, i++) {
            for (vector <Edge>::iterator it2 = dualSolution.begin(); it2 != dualSolution.end(); 
                    it2++) {
                if (areEdgesExtremesEquals((*it2), it->e) || 
                        areEdgesExtremesEquals((*it2), it->f)) {
                    G[i] += 1.0;
                }
            }

            if (fabs(u[i]) < EPSILON && G[i] < 0.0) {
                G[i] = 0.0;
            }
        }

        /* Obtaining a feasible primal solution from a (possibly unfeasible) dual solution */
        primalSolution = vector <Edge> (dualSolution);
        primalBoundValue = fixSolution(&primalSolution, n, E, S, fixedEdges, 
                fixSolutionTimeLimit);

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

        (*totalIterations)++;
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
        vector <unsigned int> * parent, unsigned int u, 
        set <Edge, bool (*) (Edge, Edge)> * bridges) {
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
                (*bridges).insert(e);
            }
        } else if (v != (*parent)[u] && (*low)[u] > (*disc)[v]) {
            (*low)[u] = (*disc)[v];
        }
    }
}

set <Edge, bool (*) (Edge, Edge)> getBridges (
        vector <list < pair <unsigned int, double> > > adj) {
    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    set <Edge, bool (*) (Edge, Edge)> result (edgeComparatorPointer);
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

bool isConflictingWithSomeEdge (
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, Edge e) {
    for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it = 
            S.begin(); it != S.end(); it++) {
        if (areEdgesExtremesEquals(e, it->e) || areEdgesExtremesEquals(e, it->f)) {
            return true;
        }
    }

    return false;
}

bool areEdgesConflicting (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        Edge e, Edge f) {
    for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it = 
            S.begin(); it != S.end(); it++) {
        if ((areEdgesExtremesEquals(e, it->e) && areEdgesExtremesEquals(f, it->f)) ||
                    (areEdgesExtremesEquals(e, it->f) && areEdgesExtremesEquals(f, it->e))) {
            return true;
        }
    }

    return false;
}

void preProcessingPhase1 (unsigned int n, vector <Edge> * E, set <ConflictingPair, 
        bool (*) (ConflictingPair, ConflictingPair)> * S, 
        set <Edge, bool (*) (Edge, Edge)> * fixedEdges, 
        chrono::high_resolution_clock::time_point tBeginPreProcessing, 
        unsigned int preProcessingTimeLimit) {
    bool flag = false;

    vector < list < pair <unsigned int, double> > > adj = toAdjacencyList (n, (*E));

    while (!timeLimitExceeded(tBeginPreProcessing, preProcessingTimeLimit) && !flag && 
            isConnected(n, (*E))) {
        set <Edge, bool (*) (Edge, Edge)> bridges = getBridges(adj);

        for (set <Edge, bool (*) (Edge, Edge)>::iterator it = (*fixedEdges).begin(); 
                it != (*fixedEdges).end(); it++) {
            bridges.erase((*it));
        }

        /* For each bridge e not yet fixed */
        for (set <Edge, bool (*) (Edge, Edge)>::iterator it = bridges.begin(); 
                it != bridges.end(); it++) {
            Edge e = (*it);

            /* Add e to the fixedEdges */
            (*fixedEdges).insert(e);

            /* For each conflicting pair containing e */
            for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it2 
                    = (*S).begin(); it2 != (*S).end();) {
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

bool preProcessingPhase2 (unsigned int n, vector <Edge> * E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> * S, 
        set <Edge, bool (*) (Edge, Edge)> * fixedEdges, 
        chrono::high_resolution_clock::time_point tBeginPreProcessing, 
        unsigned int preProcessingTimeLimit) {
    bool result = false;

    /* For each edge e that is conflicting with some edge */
    for (vector <Edge>::iterator it = (*E).begin(); 
            !timeLimitExceeded(tBeginPreProcessing, preProcessingTimeLimit) && 
            it != (*E).end();) {
        Edge e = (*it);

        if (isConflictingWithSomeEdge((*S), e)) {
            /* Create a copy of the graph */
            vector <Edge> Eprime ((*E));

            /* Create a copy of the conflicting pairs */
            set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> Sprime ((*S));

            /* Create a copy of the fixed edges */
            set <Edge, bool (*) (Edge, Edge)> fixedEdgesPrime ((*fixedEdges));

            /* Fix e in the copied instance */
            /* by removing the edges that conflict with e from the copied graph */
            /* and removing conflicts involving e from the copied conflicting pairs */
            for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it2 
                    = Sprime.begin(); it2 != Sprime.end();) {
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
            preProcessingPhase1(n, &Eprime, &Sprime, &fixedEdgesPrime, tBeginPreProcessing, 
                    preProcessingTimeLimit);

            /* If the copied graph becomes disconnected */
            if (!isConnected(n, Eprime)) {
                /* Removes conflicts involving e */
                for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator 
                        it2 = (*S).begin(); it2 != (*S).end();) {
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

bool preProcessingPhase3 (unsigned int n, vector <Edge> * E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> * S, 
        set <Edge, bool (*) (Edge, Edge)> * fixedEdges, 
        chrono::high_resolution_clock::time_point tBeginPreProcessing, 
        unsigned int preProcessingTimeLimit) {
    bool result = false;

    /* For each pair of edges (e, f) */
    /* that conflict with some edge but do not conflict with each other */
    for (vector <Edge>::iterator it = (*E).begin(); 
            !timeLimitExceeded(tBeginPreProcessing, preProcessingTimeLimit) && it != (*E).end(); 
            it++) {
        Edge e = (*it);

        if (isConflictingWithSomeEdge((*S), e)) {
            for (vector <Edge>::iterator it2 = it + 1; it2 != (*E).end(); it2++) {
                Edge f = (*it2);

                if (isConflictingWithSomeEdge((*S), f) && !areEdgesConflicting((*S), e, f)) {
                    /* Create a copy of the graph */
                    vector <Edge> Eprime ((*E));

                    /* Create a copy of the conflicting pairs */
                    set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> Sprime 
                        ((*S));

                    /* Create a copy of the fixed edges */
                    set <Edge, bool (*) (Edge, Edge)> fixedEdgesPrime ((*fixedEdges));

                    /* Fix edges e and f in the copied instance */
                    fixedEdgesPrime.insert(e);
                    fixedEdgesPrime.insert(f);

                    /* Removes from the copied instance the edges conflicting with e ou with f */
                    for (set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)>::iterator it3 = Sprime.begin(); it3 != Sprime.end();) {
                        if (areEdgesExtremesEquals(e, it3->e) ||
                                areEdgesExtremesEquals(f, it3->e)) {
                            removeEdge(&Eprime, it3->e);
                            it3 = Sprime.erase(it3);
                        } else if (areEdgesExtremesEquals(e, it3->f) ||
                                areEdgesExtremesEquals(f, it3->f)) {
                            removeEdge(&Eprime, it3->f);
                            it3 = Sprime.erase(it3);
                        } else {
                            it3++;
                        }
                    }

                    /* Apply pre-processing phase 1 in the copied instance */
                    /* in order to fix the new bridges */
                    preProcessingPhase1(n, &Eprime, &Sprime, &fixedEdgesPrime, 
                            tBeginPreProcessing, preProcessingTimeLimit);

                    /* If the copied graph becomes disconnected */
                    if (!isConnected(n, Eprime)) {
                        /* Include new conflicting pair {e, f} */
                        ConflictingPair cp;

                        if (edgeComparator(e, f)) {
                            cp.e = e;
                            cp.f = f;
                        } else {
                            cp.e = f;
                            cp.f = e;
                        }

                        (*S).insert(cp);

                        /* Update return value */
                        result = true;
                    }
                }
            }
        }
    }

    return result;
}

void preProcessing (unsigned int n, vector <Edge> * E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> * S, 
        set <Edge, bool (*) (Edge, Edge)> * fixedEdges, unsigned int preProcessingTimeLimit) {
    chrono::high_resolution_clock::time_point tBeginPreProcessing = 
        chrono::high_resolution_clock::now();

    do {
        do {
            preProcessingPhase1(n, E, S, fixedEdges, tBeginPreProcessing, preProcessingTimeLimit);
        } while (!timeLimitExceeded(tBeginPreProcessing, preProcessingTimeLimit) && 
                isConnected(n, (*E)) && 
                preProcessingPhase2(n, E, S, fixedEdges, tBeginPreProcessing, 
                    preProcessingTimeLimit));
    } while (!timeLimitExceeded(tBeginPreProcessing, preProcessingTimeLimit) && 
            isConnected(n, (*E)) && 
            preProcessingPhase3(n, E, S, fixedEdges, tBeginPreProcessing, 
                preProcessingTimeLimit));
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

    unsigned int totalTimeLimit, preProcessingTimeLimit, fixSolutionTimeLimit, N;
    double pi, minPi;

    if (!readParameters(k, &totalTimeLimit, &preProcessingTimeLimit, &fixSolutionTimeLimit, &pi, 
                &N, &minPi)) {
        cerr << "Error while reading parameters!" << endl;
        return 1;
    }

    unsigned int n;
    vector <Edge> E;
    set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S;

    if (!readInput(&n, &E, &S, argv[2])) {
        cerr << "Error while reading input!" << endl;
    }

    double bestDualBoundValue, bestPrimalBoundValue;
    int bestDualBoundIteration, totalIterations, bestPrimalBoundIteration;
    vector <Edge> bestPrimalSolution;

    bestDualBoundValue = -INFINITE;
    bestDualBoundIteration = -1;
    totalIterations = 0;
    bestPrimalBoundValue = INFINITE;
    bestPrimalBoundIteration = -1;

    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    set <Edge, bool (*) (Edge, Edge)> fixedEdges (edgeComparatorPointer);
    preProcessing (n, &E, &S, &fixedEdges, preProcessingTimeLimit);

    chrono::high_resolution_clock::time_point tCurrent = chrono::high_resolution_clock::now();
    chrono::seconds elapsedTime = chrono::duration_cast <chrono::seconds> (tCurrent - tBegin);

    if (isConnected(n, E)) {
        if (k == 1) {
            if (!relaxLag1(&bestDualBoundValue, &bestDualBoundIteration, &totalIterations, 
                        &bestPrimalBoundValue, &bestPrimalBoundIteration, &bestPrimalSolution, n, 
                        E, S, fixedEdges, tBegin, totalTimeLimit, fixSolutionTimeLimit, pi, N, 
                        minPi)) {
                cerr << "Error while executing first Lagrangian Relaxation!" << endl;
                return 1;
            }
        } else {
            if (!relaxLag2(&bestDualBoundValue, &bestDualBoundIteration, &totalIterations, 
                        &bestPrimalBoundValue, &bestPrimalBoundIteration, &bestPrimalSolution, n, 
                        E, S, fixedEdges, tBegin, totalTimeLimit, fixSolutionTimeLimit, pi, N, 
                        minPi)) {
                cerr << "Error while executing second Lagrangian Relaxation!" << endl;
                return 1;
            }
        }
    }

    if (!writeResult(argv[3], bestDualBoundValue, bestDualBoundIteration, totalIterations,
                bestPrimalBoundValue, bestPrimalBoundIteration, bestPrimalSolution)) {
        cerr << "Error while writing output!" << endl;
        return 1;
    }

    return 0;
}

