#include "graph.h"

// Graph class implementation
// Constructor 
Graph::Graph(int vertexCount) {
    this->vertexCount = vertexCount;

    // Initialize the adjacency matrix with all 0s
    adjMatrix.resize(vertexCount, std::vector<int>(vertexCount, std::numeric_limits<int>::max()));
}

void Graph::addEdge(int row, int column, int weight) {
    adjMatrix[row][column] = weight;
}

void Graph::removeEdge(int row, int column) {
    adjMatrix[row][column] = std::numeric_limits<int>::max();
}

void Graph::printMatrix() {
    for (auto row : adjMatrix) {
        for (auto column : row) {
            if (column == std::numeric_limits<int>::max())
                std::cout << "INF ";
            else
                std::cout << column << " ";
        }
        std::cout << std::endl;
    }
}

// Main function
int main() {
    Graph g(5);

    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 2);
    g.addEdge(1, 2, 3);
    g.addEdge(2, 3, 4);
    g.addEdge(3, 4, 5);

    g.printMatrix();

    return 0;
}