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
#include <map>


// Estructura para representar una arista
struct Edge {
    int src, dest, weight;
};

// Clase para representar el grafo
class Graph {
private:
    std::vector<Edge> edges;
    int vertexCount;

    // Parte 4
    std::vector<std::pair<float, float>> coordinates;

public:
    std::vector<std::vector<int>> adjMatrix, directedGraph;
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

    // 6. Calcular el área de Voronoi
    void voronoi();

    // HELPER: Remove duplicates
    void removeDuplicates(std::vector<std::pair<float, float>>& coordinates);

    std::pair<float, float> findCenter(std::vector<std::pair<float, float>> coordinates);
    bool isInside(const std::vector<std::pair<float, float>>& polygon, std::pair<float, float> point);
    bool onSegment(std::pair<float, float> p, std::pair<float, float> q, std::pair<float, float> r);
    int orientation(std::pair<float, float> p, std::pair<float, float> q, std::pair<float, float> r);
    bool doIntersect(std::pair<float, float> p1, std::pair<float, float> q1, std::pair<float, float> p2, std::pair<float, float> q2);
    
    double distanceToMidpoint(const std::pair<float, float>& point, const std::vector<std::pair<float, float>>& points);
    double distance(const std::pair<float, float>& p1, const std::pair<float, float>& p2);

private:
    int find(std::vector<int>& parent, int i);
    void unionSets(std::vector<int>& parent, std::vector<int>& rank, int x, int y);

    int bfs(std::vector<std::vector<int>>& rGraph, int s, int t, std::vector<int>& parent);
    float roundToTwoDecimals(float var);
};

#endif // GRAPH_H
