// graph.cpp

#include "graph.h"
#include <algorithm>
#include <iostream> 

// Constructor
Graph::Graph(int vertexCount) {
    this->vertexCount = vertexCount;

    adjMatrix.resize(vertexCount, std::vector<int>(vertexCount, std::numeric_limits<int>::max()));
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


int main() {
    Graph g(7);

    g.addEdge(0, 1, 3);
    g.addEdge(0, 2, 5);
    g.addEdge(0, 4, 6);
    g.addEdge(1, 3, 3);
    g.addEdge(2, 3, 2);
    g.addEdge(2, 4, 3);
    g.addEdge(2, 5, 1);
    g.addEdge(4, 5, 4);
    g.addEdge(4, 6, 2);
    g.addEdge(5, 6, 3);

    g.printMatrix();
    std::cout << std::endl;

    g.kruskalMST();

    return 0;
}

