// graph.cpp

#include "graph.h"
#include <algorithm>
#include <iostream> 

// Constructor
Graph::Graph(int vertexCount) {
    this->vertexCount = vertexCount;

    adjMatrix.resize(vertexCount, std::vector<int>(vertexCount, 0));
    directedGraph.resize(vertexCount, std::vector<int>(vertexCount, 0));
}

void Graph::addEdge(int row, int column, int weight, std::vector<std::vector<int>> &matrix) {
    matrix[row][column] = weight;
    matrix[column][row] = weight;
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

int Graph::edmondsKarp() {
    std::vector<std::vector<int>> rGraph = directedGraph;
    std::vector<int> parent(vertexCount);
    int maxFlow = 0;

    // Use the first node as the source and the last node as the sink
    int source = 0;
    int sink = vertexCount - 1;

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

        maxFlow += pathFlow;
    }

    return maxFlow;
}

/// MARK: - Voronoi, parte 4
// 1. Obtener los centroides
std::pair<float, float> Graph::getCenterPoint(std::pair<float, float> pointA, std::pair<float, float> pointB) {
    float x = (pointA.first + pointB.first) / 2;
    float y = (pointA.second + pointB.second) / 2;

    return {x, y};
}    

// 2. Obtener la pendiente
double Graph::getSlope(std::pair<float, float> pointA, std::pair<float, float> pointB) {
    return (pointB.second - pointA.second) / (pointB.first - pointA.first);
}

// 3. Calcular la mediatriz
std::pair<double, double> Graph::getMediatrix(std::pair<float, float> pointA, std::pair<float, float> pointB) {
    std::pair<float, float> center = getCenterPoint(pointA, pointB);
    double slope = getSlope(pointA, pointB);

    double perpendicularSlope = -1 / slope;
    double yIntercept = center.second - perpendicularSlope * center.first;

    return {perpendicularSlope, yIntercept};
}

// 4. Calcular la intersección
std::pair<float, float> Graph::getIntersection(std::pair<float, float> lineA, std::pair<float, float> lineB) {
    // Verificar si las líneas son paralelas
    if (lineA.first == lineB.first) {
        // Líneas paralelas, no hay intersección
        return {-1, -1}; // Otra forma de indicar que no hay intersección
    }
    float x = (lineB.second - lineA.second) / (lineA.first - lineB.first);
    float y = lineA.first * x + lineA.second;

    return {x, y};
}

// 5. Consultar el área de Voronoi
bool Graph::rayCast(const std::vector<std::pair<float, float>>& polygon, std::pair<float, float> intersection) {
    bool inside = false;
    int n = polygon.size();

    for (int i = 0, j = n - 1; i < n; j = i++) {
        if (((polygon[i].second > intersection.second) != (polygon[j].second > intersection.second)) &&
            (intersection.first < (polygon[j].first - polygon[i].first) * (intersection.second - polygon[i].second) / (polygon[j].second - polygon[i].second) + polygon[i].first)) {
            inside = !inside;
        }
    }

    return inside;
}



// 6. Calcular el área de Voronoi
void Graph::voronoi() {
    // Ordenamos los puntos en sentido horario
    for (const auto& point : coordinates) {
        std::cout << "(" << point.first << ", " << point.second << ")" << std::endl;
    }

    // Calculamos las mediatrices
    std::vector<std::pair<double, double>> mediatrixes;
    for (int i = 0; i < coordinates.size(); ++i) {
        for (int j = i + 1; j < coordinates.size(); ++j) {
            mediatrixes.push_back(getMediatrix(coordinates[i], coordinates[j]));
        }
    }

    // Calculamos las intersecciones
    std::vector<std::pair<float, float>> intersections;
    for (int i = 0; i < mediatrixes.size(); ++i) {
        for (int j = i + 1; j < mediatrixes.size(); ++j) {
            intersections.push_back(getIntersection(mediatrixes[i], mediatrixes[j]));
        }
    }

    for (const auto& intersection : intersections) {
        std::cout << "Intersección: (" << intersection.first << ", " << intersection.second << ")" << std::endl;
    }
}

/// MARK: - Sort coordinates clockwise
void Graph::sortCoordinates(std::vector<std::pair<float, float>>& coordinates) {
    float centroidX = 0, centroidY = 0;
    for (const auto& point : coordinates) {
        centroidX += point.first;
        centroidY += point.second;
    }
    centroidX /= coordinates.size();
    centroidY /= coordinates.size();

    std::sort(coordinates.begin(), coordinates.end(),
              [centroidX, centroidY](const std::pair<float, float>& a, const std::pair<float, float>& b) {
                  return atan2(a.second - centroidY, a.first - centroidX) > atan2(b.second - centroidY, b.first - centroidX);
              });
}

/// MARK: - Función para imprimir la matriz de adyacencia
void Graph::printMatrix() {
    for (auto row : adjMatrix) {
        for (auto column : row) {
            if (column == 0)
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

    std::queue<std::vector<std::vector<int>>*> matrixQueue;
    matrixQueue.push(&adjMatrix);
    matrixQueue.push(&directedGraph);

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

                for (int i = 0; i < row.size(); ++i) {
                    if (row[i] != 0 && row[i] != std::numeric_limits<int>::max()) {
                        addEdge(counter - 1, i, row[i], *matrixQueue.front());
                    }
                }

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
    g.kruskalMST();
    g.printMatrix();
    g.voronoi();

    std::cout << "Flujo máximo: " << g.edmondsKarp() << std::endl;

    return 0;
}

