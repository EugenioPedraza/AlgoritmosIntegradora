// graph.h

#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <iostream>
#include <limits>
#include <queue>
#include <utility>
#include <fstream>
#include <sstream>


// Estructura para representar una arista
struct Edge {
    int src, dest, weight;
};

// Clase para representar el grafo
class Graph {
private:
    std::vector<std::vector<int>> adjMatrix, directedGraph;
    std::vector<Edge> edges;
    int vertexCount;

    // Parte 4
    std::vector<std::pair<float, float>> coordinates;

public:
    // Constructor
    Graph(int vertexCount);

    // Métodos CRUD
    void addEdge(int row, int column, int weight, std::vector<std::vector<int>> &matrix);
    void printMatrix();

    // Encontrar MST usando Kruskal
    void kruskalMST();

    // Algoritmo de Ford Fulkerson de flujo máximo
    int edmondsKarp();

    void readFromFile(std::string filename);

    // Voronoi
    // 1. Obtener los centroides
    std::pair<float, float> getCenterPoint(std::pair<float, float> pointA, std::pair<float, float> pointB); 

    // 2. Obtener la pendiente 
    double getSlope(std::pair<float, float> pointA, std::pair<float, float> pointB);

    // 3. Calcular la mediatriz
    std::pair<double, double> getMediatrix(std::pair<float, float> pointA, std::pair<float, float> pointB);

    // 4. Calcular la intersección
    std::pair<float, float> getIntersection(std::pair<float, float> lineA, std::pair<float, float> lineB);

    // 5. Consultar el área de Voronoi
    bool rayCast(const std::vector<std::pair<float, float>>& polygon, std::pair<float, float> intersection);

    // 6. Calcular el área de Voronoi
    void voronoi();

    // HELPER: Sort coordinates clockwise
    void sortCoordinates(std::vector<std::pair<float, float>>& coordinates);

private:
    int find(std::vector<int>& parent, int i);
    void unionSets(std::vector<int>& parent, std::vector<int>& rank, int x, int y);

    int bfs(std::vector<std::vector<int>>& rGraph, int s, int t, std::vector<int>& parent);
};

#endif // GRAPH_H
