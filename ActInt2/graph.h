// graph.h

#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <iostream>
#include <limits>

// Estructura para representar una arista
struct Edge {
    int src, dest, weight;
};

// Clase para representar el grafo
class Graph {
private:
    std::vector<std::vector<int>> adjMatrix;
    std::vector<Edge> edges;
    int vertexCount;

public:
    // Constructor
    Graph(int vertexCount);

    // Métodos CRUD
    void addEdge(int row, int column, int weight);
    void printMatrix();

    // Encontrar MST usando Kruskal
    void kruskalMST();

private:
    int find(std::vector<int>& parent, int i);
    void unionSets(std::vector<int>& parent, std::vector<int>& rank, int x, int y);
};

#endif // GRAPH_H
