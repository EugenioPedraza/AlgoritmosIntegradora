// graph.cpp

#include "graph.h"
#include <algorithm>
#include <iostream> 
struct PairFloatComparator {
    bool operator()(const std::pair<float, float>& a, const std::pair<float, float>& b) const {
        float epsilon = 0.0001f; // or a suitable small number
        return std::abs(a.first - b.first) < epsilon && std::abs(a.second - b.second) < epsilon;
    }

    std::size_t operator()(const std::pair<float, float>& a) const {
        auto h1 = std::hash<float>{}(a.first);
        auto h2 = std::hash<float>{}(a.second);
        return h1 ^ h2;
    }
};

bool operator==(const std::pair<float, float>& a, const std::pair<float, float>& b) {
    float epsilon = 0.0001f; // or a suitable small number
    return std::abs(a.first - b.first) < epsilon && std::abs(a.second - b.second) < epsilon;
}

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

    x = roundToTwoDecimals(x);
    y = roundToTwoDecimals(y);

    return {x, y};
}

float Graph::roundToTwoDecimals(float var) {
    float value = (int)(var * 100 + .5);
    return (float)value / 100;
}



// 6. Calcular el área de Voronoi
void Graph::voronoi() {
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

    std::sort(coordinates.begin(), coordinates.end());
    

    // Ordenamos las intersecciones según su distancia al punto medio entre dos puntos
    std::sort(intersections.begin(), intersections.end(), [&](const auto& a, const auto& b) {
        double distA = distanceToMidpoint(a, coordinates);
        double distB = distanceToMidpoint(b, coordinates);
        return distA < distB;
    });

    std::pair<float, float> pivot = findCenter(coordinates);

    std::vector<std::pair<float, float>> filteredIntersections;
    for (const auto& intersection : intersections) {
        double distToPivot = distance(intersection, pivot);
        if (distToPivot < 10.0) { // Ajusta este valor según tu criterio
            filteredIntersections.push_back(intersection);
        }
    }

    // Almacena las intersecciones filtradas en una estructura de datos, como un mapa
    std::map<std::pair<float, float>, int> voronoiPolygons;
    for (const auto& intersection : filteredIntersections) {
        voronoiPolygons[intersection] += 1;
    }

    // Imprime las intersecciones filtradas
    for (auto element:voronoiPolygons) {
        if (element.second >= 3) {
            std::cout << "(" << element.first.first << ", " << element.first.second << ") " << element.second << std::endl;
        }
    }
}


/// MARK: - Remove duplicates
void Graph::removeDuplicates(std::vector<std::pair<float, float>>& coordinates) {
    std::sort(coordinates.begin(), coordinates.end());
    coordinates.erase(std::unique(coordinates.begin(), coordinates.end(), PairFloatComparator()), coordinates.end());
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

bool Graph::isInside(const std::vector<std::pair<float, float>>& polygon, std::pair<float, float> point) {
    int n = polygon.size();
    if (n < 3) return false;  // A polygon must have at least 3 vertices

    // Create a point for line segment from p to infinite
    std::pair<float, float> extreme = {INT_MAX, point.second};

    // Count intersections of the above line with sides of polygon
    int count = 0, i = 0;
    do {
        int next = (i + 1) % n;

        // Check if the line segment from 'p' to 'extreme' intersects
        // with the line segment from 'polygon[i]' to 'polygon[next]'
        if (doIntersect(polygon[i], polygon[next], point, extreme)) {
            // If the point 'p' is colinear with line segment 'i-next',
            // then check if it lies on segment. If it lies, return true,
            // otherwise false
            if (orientation(polygon[i], point, polygon[next]) == 0)
                return onSegment(polygon[i], point, polygon[next]);

            count++;
        }
        i = next;
    } while (i != 0);

    // Return true if count is odd, false otherwise
    return count & 1;
}

std::pair<float, float> Graph::findCenter(std::vector<std::pair<float, float>> coordinates) {
    for (int i = 0; i < coordinates.size(); ++i) {
        // Create a polygon without the i-th point
        std::vector<std::pair<float, float>> polygon;
        for (int j = 0; j < coordinates.size(); ++j) {
            if (j != i) {
                polygon.push_back(coordinates[j]);
            }
        }

        // Check if the i-th point is inside the polygon
        if (isInside(polygon, coordinates[i])) {
            return coordinates[i];
        }
    }

    // No point is inside the polygon
    return {-1, -1};
}
// Given three colinear points p, q, r, the function checks if point q lies on line segment 'pr'
bool Graph::onSegment(std::pair<float, float> p, std::pair<float, float> q, std::pair<float, float> r) {
    if (q.first <= std::max(p.first, r.first) && q.first >= std::min(p.first, r.first) &&
        q.second <= std::max(p.second, r.second) && q.second >= std::min(p.second, r.second))
        return true;
    return false;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int Graph::orientation(std::pair<float, float> p, std::pair<float, float> q, std::pair<float, float> r) {
    float val = (q.second - p.second) * (r.first - q.first) - (q.first - p.first) * (r.second - q.second);
    if (val == 0) return 0;  // colinear
    return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// The function that returns true if line segment 'p1q1' and 'p2q2' intersect.
bool Graph::doIntersect(std::pair<float, float> p1, std::pair<float, float> q1, std::pair<float, float> p2, std::pair<float, float> q2) {
    // Find the four orientations needed for general and special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and q2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

    // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases
}

double Graph::distanceToMidpoint(const std::pair<float, float>& point, const std::vector<std::pair<float, float>>& points) {
    double minDistance = std::numeric_limits<double>::max();
    for (int i = 0; i < points.size(); ++i) {
        for (int j = i + 1; j < points.size(); ++j) {
            std::pair<float, float> midpoint = getCenterPoint(points[i], points[j]);
            double dist = distance(point, midpoint);
            if (dist < minDistance) {
                minDistance = dist;
            }
        }
    }
    return minDistance;
}

double Graph::distance(const std::pair<float, float>& p1, const std::pair<float, float>& p2) {
    float dx = p2.first - p1.first;
    float dy = p2.second - p1.second;
    return sqrt(dx * dx + dy * dy);
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