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

bool isConflictingWithSomeEdge (
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, Edge e) {
    for (set <ConflictingPair>::iterator it = S.begin(); it != S.end(); it++) {
        if (areEdgesExtremesEquals(e, it->e) || areEdgesExtremesEquals(e, it->f)) {
            return true;
        }
    }

    return false;
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

vector < list < pair <unsigned int, double> > > toAdjacencyList (unsigned int n, 
        set <Edge, bool (*) (Edge, Edge)> E) {
    vector < list < pair <unsigned int, double> > > adj (n);

    for (set <Edge>::iterator it = E.begin(); it != E.end(); it++) {
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

double kruskal (set <Edge, bool (*) (Edge, Edge)> * A, unsigned int n, 
        set <Edge, bool (*) (Edge, Edge)> E) {
    double result = 0.0;

    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    (*A) = set <Edge, bool (*) (Edge, Edge)> (edgeComparatorPointer); /* A ← ∅ */

    vector <Edge> Eprime (E.begin(), E.end());

    Sets S;

    SETSinit (&S, n);

    /* For each vertex v ∈ V */
    for (unsigned int v = 0; v < n; v++) {
        SETSmake (&S, v);
    }

    /* Sort the edges of E into nondecreasing order by weight w */
    sort(Eprime.begin(), Eprime.end(), edgeWeightComparator);

    /* For each edge (u, v) ∈ E, taken in nondecreasing order by weight */
    for (vector <Edge>::iterator it = Eprime.begin(); it != Eprime.end(); it++) {
        Edge e = (*it);

        if (SETSfind (&S, e.u) != SETSfind (&S, e.v)) {
            (*A).insert(e);
            result += e.w;
            SETSunion (&S, e.u, e.v);
        }
    }

    SETSdestroy (&S);

    return result;
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

bool readParameters (unsigned int k, unsigned int * totalTimeLimit, 
        unsigned int * preProcessingTimeLimit, unsigned int * constructiveHeuristicTimeLimit, 
        unsigned int * fixSolutionTimeLimit, double * pi, unsigned int * N, double * minPi) {
    ifstream parameterFileStream(k == 1 ? "param1" : "param2", ifstream::in);

    if (!parameterFileStream.is_open()) {
        return false;
    }

    /* Total time limit in seconds */
    parameterFileStream >> (*totalTimeLimit);

    /* Pre-processing time limit in seconds */
    parameterFileStream >> (*preProcessingTimeLimit);

    /* Constructive heuristic time limit in seconds */
    parameterFileStream >> (*constructiveHeuristicTimeLimit);

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

bool readInput (unsigned int * n, set <Edge, bool (*) (Edge, Edge)> * E, set <ConflictingPair, 
        bool (*) (ConflictingPair, ConflictingPair)> * S, char * inputFilePath) {
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

bool timeLimitExceeded (chrono::high_resolution_clock::time_point tBegin, 
        unsigned int timeLimit) {
    chrono::high_resolution_clock::time_point tCurrent = chrono::high_resolution_clock::now();
    chrono::seconds elapsedTime = chrono::duration_cast <chrono::seconds> (tCurrent - tBegin);

    if ((unsigned int) elapsedTime.count() >= timeLimit) {
        return true;
    }

    return false;
}

void preProcessingPhase1 (unsigned int n, set <Edge, bool (*) (Edge, Edge)> * E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> * S, 
        set <Edge, bool (*) (Edge, Edge)> * fixedEdges, 
        chrono::high_resolution_clock::time_point tBeginPreProcessing, 
        unsigned int preProcessingTimeLimit) {
    bool flag = false;

    vector < list < pair <unsigned int, double> > > adj = toAdjacencyList (n, (*E));

    while (!timeLimitExceeded(tBeginPreProcessing, preProcessingTimeLimit) && !flag && 
            isConnected(n, (*E))) {
        set <Edge, bool (*) (Edge, Edge)> bridges = getBridges(adj);

        for (set <Edge>::iterator it = (*fixedEdges).begin(); it != (*fixedEdges).end(); it++) {
            bridges.erase((*it));
        }

        /* For each bridge e not yet fixed */
        for (set <Edge>::iterator it = bridges.begin(); it != bridges.end(); it++) {
            Edge e = (*it);

            /* Add e to the fixedEdges */
            (*fixedEdges).insert(e);

            /* For each conflicting pair containing e */
            for (set <ConflictingPair>::iterator it2 = (*S).begin(); it2 != (*S).end();) {
                ConflictingPair cp = (*it2);

                if (areEdgesExtremesEquals(e, cp.e)) {
                    /* Remove the conflicing edge from the graph */
                    (*E).erase(cp.f);
                    removeEdgeFromAdjacencyList(adj, cp.f);

                    /* Remove the conflicting pair */
                    it2 = (*S).erase(it2);
                } else if (areEdgesExtremesEquals(e, cp.f)) {
                    /* Remove the conflicing edge from the graph */
                    (*E).erase(cp.e);
                    removeEdgeFromAdjacencyList(adj, cp.e);

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

bool preProcessingPhase2 (unsigned int n, set <Edge, bool (*) (Edge, Edge)> * E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> * S, 
        set <Edge, bool (*) (Edge, Edge)> * fixedEdges, 
        chrono::high_resolution_clock::time_point tBeginPreProcessing, 
        unsigned int preProcessingTimeLimit) {
    bool result = false;

    /* For each edge e that is conflicting with some edge */
    for (set <Edge>::iterator it = (*E).begin(); 
            !timeLimitExceeded(tBeginPreProcessing, preProcessingTimeLimit) && 
            it != (*E).end();) {
        Edge e = (*it);

        if (isConflictingWithSomeEdge((*S), e)) {
            /* Create a copy of the graph */
            set <Edge, bool (*) (Edge, Edge)> Eprime ((*E));

            /* Create a copy of the conflicting pairs */
            set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> Sprime ((*S));

            /* Create a copy of the fixed edges */
            set <Edge, bool (*) (Edge, Edge)> fixedEdgesPrime ((*fixedEdges));

            /* Fix e in the copied instance */
            fixedEdgesPrime.insert(e);

            /* Remove the edges that conflicts with e from the copied graph */
            /* and remove conflicts involving e from the copied conflicting pairs */
            for (set <ConflictingPair>::iterator it2 = Sprime.begin(); it2 != Sprime.end();) {
                ConflictingPair cp = (*it2);

                if (areEdgesExtremesEquals(e, cp.e)) {
                    Eprime.erase(cp.f);

                    it2 = Sprime.erase(it2);
                } else if (areEdgesExtremesEquals(e, cp.f)) {
                    Eprime.erase(cp.e);

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
                for (set <ConflictingPair>::iterator it2 = (*S).begin(); it2 != (*S).end();) {
                    ConflictingPair cp = (*it2);

                    if (areEdgesExtremesEquals(e, cp.e) || areEdgesExtremesEquals(e, cp.f)) {
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

bool preProcessingPhase3 (unsigned int n, set <Edge, bool (*) (Edge, Edge)> * E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> * S, 
        set <Edge, bool (*) (Edge, Edge)> * fixedEdges, 
        chrono::high_resolution_clock::time_point tBeginPreProcessing, 
        unsigned int preProcessingTimeLimit) {
    bool result = false;

    /* For each pair of edges (e, f) */
    /* that conflict with some edge but do not conflict with each other */
    for (set <Edge>::iterator it = (*E).begin(); 
            !timeLimitExceeded(tBeginPreProcessing, preProcessingTimeLimit) && it != (*E).end(); 
            it++) {
        Edge e = (*it);

        if (isConflictingWithSomeEdge((*S), e)) {
            for (set <Edge>::iterator it2 = next(it); it2 != (*E).end(); it2++) {
                Edge f = (*it2);

                ConflictingPair cp;

                if (edgeComparator(e, f)) {
                    cp.e = e;
                    cp.f = f;
                } else {
                    cp.e = f;
                    cp.f = e;
                }

                if (isConflictingWithSomeEdge((*S), f) && (*S).find(cp) == (*S).end()) {
                    /* Create a copy of the graph */
                    set <Edge, bool (*) (Edge, Edge)> Eprime ((*E));

                    /* Create a copy of the conflicting pairs */
                    set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> Sprime ((
                                *S));

                    /* Create a copy of the fixed edges */
                    set <Edge, bool (*) (Edge, Edge)> fixedEdgesPrime ((*fixedEdges));

                    /* Fix edges e and f in the copied instance */
                    fixedEdgesPrime.insert(e);
                    fixedEdgesPrime.insert(f);

                    /* Removes from the copied instance the edges conflicting with e ou with f */
                    /* and remove conflicts involving e or f from the copied conflicting pairs */
                    for (set <ConflictingPair>::iterator it3 = Sprime.begin(); 
                            it3 != Sprime.end();) {
                        ConflictingPair cp = (*it3);

                        if (areEdgesExtremesEquals(e, cp.e) || areEdgesExtremesEquals(f, cp.e)) {
                            Eprime.erase(cp.e);
                            it3 = Sprime.erase(it3);
                        } else if (areEdgesExtremesEquals(e, cp.f) || 
                                areEdgesExtremesEquals(f, cp.f)) {
                            Eprime.erase(cp.f);
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

void preProcessing (unsigned int n, set <Edge, bool (*) (Edge, Edge)> * E, 
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

double fixSolution (set <Edge, bool (*) (Edge, Edge)> * primalSolution, unsigned int n, 
        set <Edge, bool (*) (Edge, Edge)> E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        set <Edge, bool (*) (Edge, Edge)> fixedEdges, unsigned int fixSolutionTimeLimit) {
    double primalBoundValue = 0.0;

    chrono::high_resolution_clock::time_point tBeginFixSolution = 
        chrono::high_resolution_clock::now();

    /* Fixing edges costs */
    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    set <Edge, bool (*) (Edge, Edge)> primalSolutionPrime (edgeComparatorPointer);

    for (set <Edge>::iterator it = (*primalSolution).begin(); it != (*primalSolution).end() && 
            !timeLimitExceeded(tBeginFixSolution, fixSolutionTimeLimit); it++) {
        Edge e = (*it);

        set <Edge, bool (*) (Edge, Edge)>::iterator it2 = E.find(e);

        if (it2 != E.end()) {
            Edge f = (*it2);

            primalSolutionPrime.insert(f);

            primalBoundValue += f.w;
        }
    }

    (*primalSolution) = set <Edge, bool (*) (Edge, Edge)> (primalSolutionPrime);

    /* Compute edges in primal solution that are conflicting and are not fixed */
    set <Edge, bool (*) (Edge, Edge)> edgesToRemove (edgeComparatorPointer);

    for (set <ConflictingPair>::iterator it = S.begin(); it != S.end() && 
            !timeLimitExceeded(tBeginFixSolution, fixSolutionTimeLimit); it++) {
        ConflictingPair cp = (*it);

        if ((*primalSolution).find(cp.e) != (*primalSolution).end() && 
                (*primalSolution).find(cp.f) != (*primalSolution).end()) {
            if (fixedEdges.find(cp.e) == fixedEdges.end()) {
                edgesToRemove.insert(cp.e);
            }

            if (fixedEdges.find(cp.f) == fixedEdges.end()) {
                edgesToRemove.insert(cp.f);
            }
        }
    }

    /* Remove non-fixed conflicting edges from primal solution */
    for (set <Edge>::iterator it = edgesToRemove.begin(); it != edgesToRemove.end(); it++) {
        Edge e = (*it);

        (*primalSolution).erase(e);
        primalBoundValue -= e.w;
    }

    /* Compute the eges in E that are not in primal solution */
    /* nor conflicts with any edge in the primal solution */
    vector <Edge> edgesToInsert;

    for (set <Edge>::iterator it = E.begin(); it != E.end() && 
            !timeLimitExceeded(tBeginFixSolution, fixSolutionTimeLimit); it++) {
        Edge e = (*it);

        if ((*primalSolution).find(e) == (*primalSolution).end()) {
            bool isConflicting = false;

            for (set <Edge>::iterator it2 = (*primalSolution).begin(); 
                    it2 != (*primalSolution).end() && !isConflicting && 
                    !timeLimitExceeded(tBeginFixSolution, fixSolutionTimeLimit); it2++) {
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
                    isConflicting = true;
                }
            }

            if (!isConflicting) {
                edgesToInsert.push_back(e);
            }
        }
    }

    Sets components = connectedComponents(n, (*primalSolution));

    /* Sort the edges to insert into nondecreasing order by weight w */
    sort(edgesToInsert.begin(), edgesToInsert.end(), edgeWeightComparator);

    /* For each edge to insert, taken in nondecreasing order by weight */
    for (vector <Edge>::iterator it = edgesToInsert.begin(); it != edgesToInsert.end() && 
            !timeLimitExceeded(tBeginFixSolution, fixSolutionTimeLimit); it++) {
        Edge e = (*it);

        /* Check if e conflicts if some edge in the primal solution */
        bool conflicts = false;
        for (set <Edge>::iterator it2 = (*primalSolution).begin(); 
                it2 != (*primalSolution).end() && !conflicts; it2++) {
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
                conflicts = true;
            }
        }

        /* If e does not conflict with the primal solution and */
        /* connects two distinct connected components of the primal solution */
        if (!conflicts && !sameComponent(components, e.u, e.v)) {
            /* Insert e into the primal solution */
            (*primalSolution).insert(e);
            primalBoundValue += e.w;

            /* Connect the connected components */
            SETSunion(&components, e.u, e.v);
        }
    }

    SETSdestroy(&components);

    return primalBoundValue;
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

double constructiveHeuristic (set <Edge, bool (*) (Edge, Edge)> * primalSolution, unsigned int n, 
        set <Edge, bool (*) (Edge, Edge)> E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        unsigned int constructiveHeuristicTimeLimit) {
    chrono::high_resolution_clock::time_point tBeginConstructiveHeuristic = 
        chrono::high_resolution_clock::now();

    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    set <Edge, bool (*) (Edge, Edge)> conflictingEdges (edgeComparatorPointer);

    for (set <ConflictingPair>::iterator it = S.begin(); it != S.end(); it++) {
        ConflictingPair cp = (*it);

        conflictingEdges.insert(cp.e);
        conflictingEdges.insert(cp.f);
    }

    /* While S = ∅ */
    while (S.size() > 0 && 
            !timeLimitExceeded(tBeginConstructiveHeuristic, constructiveHeuristicTimeLimit)) {
        unsigned int maxConflictsCounter = 0;
        Edge maxConflictingEdge;

        /* In the graph G we choose an e ∈ E */
        /* which appears in the maximum number of conflict pairs */
        for (set <Edge>::iterator it = conflictingEdges.begin(); it != conflictingEdges.end() && 
                !timeLimitExceeded(tBeginConstructiveHeuristic, constructiveHeuristicTimeLimit); 
                it++) {
            Edge e = (*it);
            unsigned int conflictsCounter = 0;

            for (set <ConflictingPair>::iterator it2 = S.begin(); 
                    it2 != S.end() && 
                    !timeLimitExceeded(tBeginConstructiveHeuristic, 
                        constructiveHeuristicTimeLimit); 
                    it2++) {
                ConflictingPair cp = (*it2);

                if (areEdgesExtremesEquals(e, cp.e) || areEdgesExtremesEquals(e, cp.f)) {
                    conflictsCounter++;
                }
            }

            if (maxConflictsCounter < conflictsCounter) {
                maxConflictsCounter = conflictsCounter;
                maxConflictingEdge = (*it);
            }
        }

        /* Delete e from G */
        E.erase(maxConflictingEdge);

        /* And delete from S all the conflict pairs containing e */
        for (set <ConflictingPair>::iterator it = S.begin(); it != S.end() && 
                !timeLimitExceeded(tBeginConstructiveHeuristic, 
                    constructiveHeuristicTimeLimit);) {
            ConflictingPair cp = (*it);

            if (areEdgesExtremesEquals(maxConflictingEdge, cp.e) || 
                    areEdgesExtremesEquals(maxConflictingEdge, cp.f)) {
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

vector <double> initialLagrangeMultipliers (int m) {
    vector <double> result (m, 1.0);
    return result;
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

bool relaxLag1 (double * bestDualBoundValue, int * bestDualBoundIteration, int * totalIterations, 
        double * bestPrimalBoundValue, int * bestPrimalBoundIteration, 
        set <Edge, bool (*) (Edge, Edge)> * bestPrimalSolution, unsigned int n, 
        set <Edge, bool (*) (Edge, Edge)> E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        set <Edge, bool (*) (Edge, Edge)> fixedEdges, 
        chrono::high_resolution_clock::time_point tBegin, unsigned int totalTimeLimit, 
        unsigned int constructiveHeuristicTimeLimit, unsigned int fixSolutionTimeLimit, 
        double pi, unsigned int N, double minPi) {
    vector <double> u;
    unsigned int iterationsWithoutImprovment;
    double dualBoundValue, primalBoundValue;

    set <Edge, bool (*) (Edge, Edge)> dualSolution, primalSolution;

    dualBoundValue = kruskal(&dualSolution, n, E);

    if ((*bestDualBoundValue) < dualBoundValue) {
        (*bestDualBoundValue) = dualBoundValue;
    }

    primalSolution = set <Edge, bool (*) (Edge, Edge)> (dualSolution);

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

    constructiveHeuristic(&primalSolution, n, E, S, constructiveHeuristicTimeLimit);

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
        bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
        set <Edge, bool (*) (Edge, Edge)> Eu(edgeComparatorPointer);
        vector <double> G(S.size(), -1.0);

        /* Solving the Lagrangian problem with the current set of multipliers */
        for (set <Edge>::iterator it = E.begin(); it != E.end(); it++) {
            Edge e = (*it);

            unsigned int i = 0;
            for (set <ConflictingPair>::iterator it2 = S.begin(); it2 != S.end(); it2++, i++) {
                if (areEdgesExtremesEquals(e, it2->e) || areEdgesExtremesEquals(e, it2->f)) {
                    e.w += u[i];
                }
            }

            Eu.insert(e);
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
        unsigned i = 0;
        for (set <ConflictingPair>::iterator it = S.begin(); it != S.end(); it++, i++) {
            if (dualSolution.find(it->e) != dualSolution.end()) {
                G[i] += 1.0;
            }

            if (dualSolution.find(it->f) != dualSolution.end()) {
                G[i] += 1.0;
            }

            if (fabs(u[i]) < EPSILON && G[i] < 0.0) {
                G[i] = 0.0;
            }
        }

        /* Obtaining a feasible primal solution from a (possibly unfeasible) dual solution */
        primalSolution = set <Edge, bool (*) (Edge, Edge)> (dualSolution);
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
        set <Edge, bool (*) (Edge, Edge)> * bestPrimalSolution, unsigned int n, 
        set <Edge, bool (*) (Edge, Edge)> E, 
        set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S, 
        set <Edge, bool (*) (Edge, Edge)> fixedEdges, 
        chrono::high_resolution_clock::time_point tBegin, unsigned int totalTimeLimit, 
        unsigned int constructiveHeuristicTimeLimit, unsigned int fixSolutionTimeLimit, 
        double pi, unsigned int N, double minPi) {
    vector < set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> > Se;
    unsigned int eStarIndex;
    Edge eStar;
    set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> SminusSeStar (S);
    bool (*edgeComparatorPointer) (Edge, Edge) = edgeComparator;
    set <Edge, bool (*) (Edge, Edge)> EeStar (edgeComparatorPointer), 
        EminusEeStar (edgeComparatorPointer), dualSolution, primalSolution;
    vector <double> u;
    unsigned int iterationsWithoutImprovment;
    double dualBoundValue, primalBoundValue;

    dualBoundValue = kruskal(&dualSolution, n, E);

    if ((*bestDualBoundValue) < dualBoundValue) {
        (*bestDualBoundValue) = dualBoundValue;
    }

    primalSolution = set <Edge, bool (*) (Edge, Edge)> (dualSolution);
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

    constructiveHeuristic(&primalSolution, n, E, S, constructiveHeuristicTimeLimit);

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
    bool (*conflictingPairComparatorPointer) (ConflictingPair, ConflictingPair) = 
        conflictingPairComparator;
    Se = vector < set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> > (E.size(), 
            set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> (
                conflictingPairComparatorPointer));

    unsigned int i = 0;
    for (set <Edge>::iterator it = E.begin(); it != E.end(); it++, i++) {
        Edge e = (*it);

        for (set <ConflictingPair>::iterator it2 = S.begin(); it2 != S.end(); it2++) {
            ConflictingPair cp = (*it2);

            if (areEdgesExtremesEquals(e, cp.e) || areEdgesExtremesEquals(e, cp.f)) {
                Se[i].insert(cp);
            }
        }
    }

    /* Seja ainda e* a aresta de E com menor custo para a qual Se não é vazio */
    eStarIndex = i = 0;
    eStar.w = INFINITE;
    for (set <Edge>::iterator it = E.begin(); it != E.end(); it++, i++) {
        Edge e = (*it);

        if (Se[i].size() > 0 && eStar.w > e.w) {
            eStarIndex = i;
            eStar = e;
        }
    }

    for (set <ConflictingPair>::iterator it = Se[eStarIndex].begin(); it != Se[eStarIndex].end(); 
            it++) {
        ConflictingPair cp = (*it);

        SminusSeStar.erase(cp);
    }

    for (set <Edge>::iterator it = E.begin(); it != E.end(); it++) {
        Edge e = (*it);

        if (!areEdgesExtremesEquals(e, eStar)) {
            ConflictingPair cp;

            if (edgeComparator(e, eStar)) {
                cp.e = e;
                cp.f = eStar;
            } else {
                cp.e = eStar;
                cp.f = e;
            }

            if (Se[eStarIndex].find(cp) != Se[eStarIndex].end()) {
                EeStar.insert(e);
            }
        }
    }

    for (set <Edge>::iterator it = E.begin(); it != E.end(); it++) {
        Edge e = (*it);

        if (EeStar.find(e) == EeStar.end()) {
            EminusEeStar.insert((*it));
        }
    }

    u = initialLagrangeMultipliers(SminusSeStar.size());
    iterationsWithoutImprovment = 0;

    while (!termination(tBegin, totalTimeLimit, (*bestDualBoundValue), (*bestPrimalBoundValue), 
                pi, minPi)) {
        double stepSize;
        set <Edge, bool (*) (Edge, Edge)> Eu(edgeComparatorPointer);
        vector <double> G(SminusSeStar.size(), -1.0);

        /* Solving the Lagrangian problem with the current set of multipliers */
        for (set <Edge>::iterator it = EminusEeStar.begin(); it != EminusEeStar.end(); it++) {
            Edge e = (*it);

            unsigned int i = 0;
            for (set <ConflictingPair>::iterator it2 = SminusSeStar.begin(); 
                    it2 != SminusSeStar.end(); it2++, i++) {
                if (areEdgesExtremesEquals(e, it2->e) || areEdgesExtremesEquals(e, it2->f)) {
                    e.w += u[i];
                }
            }

            Eu.insert(e);
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
        unsigned int i = 0;
        for (set <ConflictingPair>::iterator it = SminusSeStar.begin(); it != SminusSeStar.end(); 
                it++, i++) {
            ConflictingPair cp = (*it);

            if (dualSolution.find(cp.e) != dualSolution.end()) {
                G[i] += 1.0;
            }

            if (dualSolution.find(cp.f) != dualSolution.end()) {
                G[i] += 1.0;
            }

            if (fabs(u[i]) < EPSILON && G[i] < 0.0) {
                G[i] = 0.0;
            }
        }

        /* Obtaining a feasible primal solution from a (possibly unfeasible) dual solution */
        primalSolution = set <Edge, bool (*) (Edge, Edge)> (dualSolution);
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
        set <Edge, bool (*) (Edge, Edge)> bestPrimalSolution) {
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
    for (set <Edge>::iterator it = bestPrimalSolution.begin(); it != bestPrimalSolution.end(); 
            it++) {
        Edge e = (*it);

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

    unsigned int totalTimeLimit, preProcessingTimeLimit, constructiveHeuristicTimeLimit, 
                 fixSolutionTimeLimit, N;
    double pi, minPi;

    if (!readParameters(k, &totalTimeLimit, &preProcessingTimeLimit, 
                &constructiveHeuristicTimeLimit, &fixSolutionTimeLimit, &pi, &N, &minPi)) {
        cerr << "Error while reading parameters!" << endl;
        return 1;
    }

    unsigned int n;
    set <Edge, bool (*) (Edge, Edge)> E;
    set <ConflictingPair, bool (*) (ConflictingPair, ConflictingPair)> S;

    if (!readInput(&n, &E, &S, argv[2])) {
        cerr << "Error while reading input!" << endl;
    }

    double bestDualBoundValue, bestPrimalBoundValue;
    int bestDualBoundIteration, totalIterations, bestPrimalBoundIteration;
    set <Edge, bool (*) (Edge, Edge)> bestPrimalSolution;

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
                        E, S, fixedEdges, tBegin, totalTimeLimit, constructiveHeuristicTimeLimit, 
                        fixSolutionTimeLimit, pi, N, minPi)) {
                cerr << "Error while executing first Lagrangian Relaxation!" << endl;
                return 1;
            }
        } else {
            if (!relaxLag2(&bestDualBoundValue, &bestDualBoundIteration, &totalIterations, 
                        &bestPrimalBoundValue, &bestPrimalBoundIteration, &bestPrimalSolution, n, 
                        E, S, fixedEdges, tBegin, totalTimeLimit, constructiveHeuristicTimeLimit, 
                        fixSolutionTimeLimit, pi, N, minPi)) {
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

