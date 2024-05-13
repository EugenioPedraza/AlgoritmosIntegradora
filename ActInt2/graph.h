#include <iostream>
#include <vector>
#include <limits>


class Graph {
    private:
    std::vector<std::vector<int>> adjMatrix;
    int vertexCount;

    public:
    // Constructor
    Graph(int vertexCount);

    // CRUD methods
    void addEdge(int row, int column, int weight);

    void removeEdge(int row, int column);

    // Print the adjacency matrix
    void printMatrix();

};
