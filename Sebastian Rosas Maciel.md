# Sebastián Rosas Maciel - A01233989
## Reflexión: Actividad Integradora No. 1

### Algoritmo de Manacher - Complejidad: O(n)
Mi mayor contribución dentro de esta actividad fue la segunda parte del reto, la cual consiste en buscar palindromos dentro de el archivo, en un principio considere que la mejor manera de encontrar codigo espejeado dentro del string era mediante recorridos dentro de un suffix tree o un suffix array, al final termine optando por utilizar el algoritmo de Manacher ya que este busca encontrar el palindromo mas largo posible.

```cpp
std::pair<int, int> Algorithm::manacher(std::string transmissionText){ 
    int stringSize = transmissionText.size();

    if (stringSize == 0)
        return {-1, -1};

    stringSize = 2 * stringSize + 1; // Position count
    std::vector<int> lpsLengthArray(stringSize, 0); // LPS Length Array
    lpsLengthArray[0] = 0;
    lpsLengthArray[1] = 1;
    int centerPosition = 1; // centerPosition
    int centerRightPosition = 2; // centerRightPosition
    int i = 0; // currentRightPosition
    int iMirror; // currentLeftPosition
    int maxLPSLength = 0;
    int maxLPSCenterPosition = 0;
    int start = -1;
    int end = -1;
    int diff = -1;

    for (i = 2; i < stringSize; i++) {
        // get currentLeftPosition iMirror for currentRightPosition i
        iMirror = 2*centerPosition-i;
        lpsLengthArray[i] = 0;
        diff = centerRightPosition - i;
        // If currentRightPosition i is within centerRightPosition R
        if (diff > 0)
            lpsLengthArray[i] = std::min(lpsLengthArray[iMirror], diff);

        // Expand palindrome centered at currentRightPosition i
        while (((i + lpsLengthArray[i]) < stringSize && (i - lpsLengthArray[i]) > 0) && 
            (((i + lpsLengthArray[i] + 1) % 2 == 0) || 
            (transmissionText[(i + lpsLengthArray[i] + 1)/2] == transmissionText[(i - lpsLengthArray[i] - 1)/2]))) {
            lpsLengthArray[i]++;
        }

        if (lpsLengthArray[i] > maxLPSLength) { // Track maxLPSLength
            maxLPSLength = lpsLengthArray[i];
            maxLPSCenterPosition = i;
        }

        // If palindrome centered at currentRightPosition i 
        // expand centerRightPosition R,
        // adjust centerPosition C based on expanded palindrome.
        if (i + lpsLengthArray[i] > centerRightPosition) {
            centerPosition = i;
            centerRightPosition = i + lpsLengthArray[i];
        }
    }
    start = (maxLPSCenterPosition - maxLPSLength)/2;
    end = start + maxLPSLength - 1;
    return {start/2, end/2}; // Return the start and end index in the original string
}
```

### Funcionamiento 
El algoritmo toma a un solo string como entrada, y regresa a un par de **dos numeros enteros** como salida, donde el primer numero representa al **indice donde empieza** el palindromo, y el segundo numero marca al indice donde **finaliza**.

Empezaremos por obtener la longitud del string introducido, comparamos su longitud con **0** para determinar si el string esta vacio, en caso de que el string este vacio, se arroja un par **{-1, -1}**, simbolizando que el string esta vacio.
```cpp
std::pair<int, int> Algorithm::manacher(std::string transmissionText){ 
    int N = transmissionText.size();

    if (N == 0)
        return {-1, -1};
```
Duplicamos el tamaño del string y le sumamos una unidad, esta modificacion sirve para manejar palindromos de longitud par.
Creamos un vector que se ajuste al nuevo tamaño.
```cpp
stringSize = 2 * stringSize + 1;
Position count
    std::vector<int> lpsLengthArray(stringSize, 0); // LPS Length Array
```