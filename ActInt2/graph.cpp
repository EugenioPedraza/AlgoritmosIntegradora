// graph.cpp

#include "graph.h"
#include <algorithm>
#include <iostream> 

// Constructor
Graph::Graph(int vertexCount) {
    this->vertexCount = vertexCount;

    adjMatrix.resize(vertexCount, std::vector<int>(vertexCount, std::numeric_limits<int>::max()));
    directedGraph.resize(vertexCount, std::vector<int>(vertexCount, std::numeric_limits<int>::max()));
}

void Graph::addEdge(int row, int column, int weight) {
    adjMatrix[row][column] = weight;
    adjMatrix[column][row] = weight;
    edges.push_back({row, column, weight});
}

// Función para encontrar el conjunto de un elemento i 
int Graph::find(std::vector<int>& parent, int i) {
    if (parent[i] != i)
        parent[i] = find(parent, parent[i]);
    return parent[i];
}

// Función para unir dos conjuntos x e y 
void Graph::unionSets(std::vector<int>& parent, std::vector<int>& rank, int x, int y) {
    int rootX = find(parent, x);
    int rootY = find(parent, y);

    if (rank[rootX] < rank[rootY])
        parent[rootX] = rootY;
    else if (rank[rootX] > rank[rootY])
        parent[rootY] = rootX;
    else {
        parent[rootY] = rootX;
        rank[rootX]++;
    }
}

/// MARK: - Kruskal, parte 1
// Algoritmo de Kruskal para encontrar el MST
void Graph::kruskalMST() {
    std::vector<Edge> result;
    int mstWeight = 0;

    // Ordenar las aristas por peso
    std::sort(edges.begin(), edges.end(), [](Edge a, Edge b) {
        return a.weight < b.weight;
    });

    std::vector<int> parent(vertexCount);
    std::vector<int> rank(vertexCount, 0);

    // Inicializar los conjuntos disjuntos
    for (int v = 0; v < vertexCount; ++v)
        parent[v] = v;

    for (const auto& edge : edges) {
        int u = find(parent, edge.src);
        int v = find(parent, edge.dest);

        if (u != v) {
            result.push_back(edge);
            mstWeight += edge.weight;
            unionSets(parent, rank, u, v);
        }
    }

    std::cout << "Costo=" << mstWeight << std::endl;
    for (const auto& edge : result) {
        std::cout << char('A' + edge.src) << "-" << char('A' + edge.dest) << std::endl;
    }
}

/// MARK: - Algoritmo de Edmonds Karp, parte 3
int Graph::bfs(std::vector<std::vector<int>>& rGraph, int s, int t, std::vector<int>& parent) {
    // Checlist to check if a node has already been visited, all fields are initialized to false
    std::vector<bool> visited(rGraph.size(), false);
    std::queue<int> q;

    q.push(s);
    visited[s] = true;
    parent[s] = -1;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int v = 0; v < vertexCount; ++v) {
            if (!visited[v] && rGraph[u][v] > 0) {
                if (v == t) {
                    parent[v] = u;
                    return true;
                }

                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }

    return false;
}

int Graph::edmondsKarp(int source, int sink) {
    std::vector<std::vector<int>> rGraph = adjMatrix;
    std::vector<int> parent(vertexCount);
    int maxFlow = 0;

    while (bfs(rGraph, source, sink, parent)) {
        int pathFlow = std::numeric_limits<int>::max();

        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            pathFlow = std::min(pathFlow, rGraph[u][v]);
        }

        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            rGraph[u][v] -= pathFlow;
            rGraph[v][u] += pathFlow;
        }
    }

    return maxFlow;
}

/// MARK: - Función para imprimir la matriz de adyacencia
void Graph::printMatrix() {
    for (auto row : adjMatrix) {
        for (auto column : row) {
            if (column == std::numeric_limits<int>::max())
                std::cout << "0 ";
            else
                std::cout << column << " ";
        }

        std::cout << std::endl;
    }
}

/// MARK: - Read from file function
void Graph::readFromFile(std::string filename) {
    std::cout << "Leyendo archivo " << filename << std::endl;
    std::ifstream input(filename);
    std::string line;

    std::queue<std::vector<std::vector<int>>> matrixQueue;
    matrixQueue.push(adjMatrix);
    matrixQueue.push(directedGraph);

    if (input.is_open()) {
        int counter = 0;
        while (std::getline(input, line)) {
            counter++;

            if (line.empty()) {
                counter = 0;
                matrixQueue.pop();

                std::cout << std::endl;

                continue;
            }

            if (matrixQueue.empty()) {
                std::istringstream iss(line);
                int first, second;
                char extraCharacters;

                if (iss >> extraCharacters >> first >> extraCharacters >> second >> extraCharacters) {
                    std::cout << "(" << first << "," << second << ")" << std::endl;
                    coordinates.push_back({first, second});
                } else {
                    std::cerr << "Error al leer el archivo" << std::endl;
                }

            } else {
                std::istringstream iss(line);
                std::vector<int> row;
                int value;

                while (iss >> value) {
                    std::cout << value << " ";
                    row.push_back(value);
                }

                matrixQueue.front()[counter - 1] = row;
                std::cout << std::endl;
            }
        }
    } else {
        std::cerr << "Error al abrir el archivo" << std::endl;
    }
}


int main() {
    Graph g(7);

    g.readFromFile("in.txt");

    return 0;
}

