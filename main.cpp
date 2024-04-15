#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "res.h"

// Part 1, find the Mcode file inside transmission text
std::vector<int> Algorithm::computeKMPTable(const std::string& word) {
    int m = word.length();
    std::vector<int> kmpTable(m, 0);
    int len = 0;
    int i = 1;
    
    while (i < m) {
        if (word[i] == word[len]) {
            len++;
            kmpTable[i] = len;
            i++;
        } else {
            if (len != 0) {
                len = kmpTable[len - 1];
            } else {
                   kmpTable[i] = 0;
                    i++;
            }
        }
    }

    return kmpTable;
}

std::tuple<bool,int, int> Algorithm::KMP_Search(const std::string& text, const std::string& pattern) {
    std::vector<int> kmpTable = computeKMPTable(pattern);
    int i = 0;
    int j = 0; 

    while (i < text.size()) {
        if (pattern[j] == text[i]) {
            j++;
            i++;
        }

        if (j == pattern.size()) {

            return std::make_tuple(true, i-j, i-1);

            j = kmpTable[j - 1];
        } else if (i < text.size() && pattern[j] != text[i]) {
            if (j != 0) {
                j = kmpTable[j - 1];
            } else {
                i = i + 1;
            }
        }
    }

    return std::make_tuple(false, -1, -1); 
}


    
std::string Algorithm::readFileIntoString(const std::string& path) {
    std::ifstream inputFile(path);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening file " << path << std::endl;
        return "";
    }
    std::stringstream buffer;
    buffer << inputFile.rdbuf();
    inputFile.close();
    return buffer.str();
}

// Part 2, look for palindromes 
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
        while (((i + lpsLengthArray[i]) < stringSize && (i - lpsLengthArray[i]) > 0)) {
            int rightIndex = (i + lpsLengthArray[i] + 1)/2;
            int leftIndex = (i - lpsLengthArray[i] - 1)/2;

            if (rightIndex >= transmissionText.size() || leftIndex < 0) {
                break;
            }

            if ((rightIndex + leftIndex) % 2 == 0 || 
                transmissionText[rightIndex] == transmissionText[leftIndex]) {
                lpsLengthArray[i]++;
            } else {
                break;
            }
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

// Part 3, find coincidences vbetween transmission files
std::vector<int> Algorithm::computeLPS(std::string pattern) {
    int M = pattern.length();
    std::vector<int> lps(M);
    lps[0] = 0;
    int len = 0;
    int i = 1;

    while (i < M) {
        if (pattern[i] == pattern[len]) {
            len++;
            lps[i] = len;
            i++;
        }
        else {
            if (len != 0) {
                len = lps[len - 1];
            }
            else {
                lps[i] = 0;
                i++;
            }
        }
    }

    return lps;
}

void Algorithm::lcs(std::string trans1, std::string trans2)
{
    int n = trans1.length();
    int m = trans2.length();
    std::vector<std::vector<int>> dp(n + 1,std::vector<int>(m + 1, 0));
    int maxLenght = 0;
    int endIndexTrans1 = 0;

    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < m; j++)
        {
            if (trans1[i-1] == trans2[j-1])
            {
                dp[i][j] = dp[i - 1][j - 1] + 1;
                if (dp[i][j]>maxLenght)
                {
                    maxLenght = dp[i][j];
                    endIndexTrans1 = i;

                }
            }
        }
    }
    if (maxLenght==0)
    {
        std::cout << "No hay coincidencia" << std::endl;
        return;
    }
    else
    {
        int starIndexTrans1 = (endIndexTrans1 - maxLenght) + 1;
        std::cout << "Se encontro coincidencia desde " << starIndexTrans1 << " hasta " << endIndexTrans1 << std::endl;
    }

}

int main() {
    Algorithm buscador;

    std::string transmission1 = buscador.readFileIntoString("transmission1V2.txt");
    std::string transmission2 = buscador.readFileIntoString("transmission2V2.txt");
    std::string mcode1 = buscador.readFileIntoString("mcode1V2.txt");

    std::cout << "\nParte 1: Busqueda de Mcode dentro de los ashivos de transmicion" << std::endl;
    std::tuple<bool,int, int> tr1_mc1 = buscador.KMP_Search(transmission1, mcode1);
    std::tuple<bool, int, int> tr2_mc1 = buscador.KMP_Search(transmission2, mcode1);

    if (std::get<0>(tr1_mc1)) {
        std::cout << "Mcode1 se encuentra dentro de Transmision1 en la posicion: " 
                  << std::get<1>(tr1_mc1) << " a " << std::get<2>(tr1_mc1) << std::endl;
    } else {
        std::cout << "Mcode1 no se encuentra dentro de Transmision1" << std::endl;
    }

    if (std::get<0>(tr2_mc1)) {
        std::cout << "Mcode1 se encuentra dentro de Transmision2 en la posicion: " 
                  << std::get<1>(tr2_mc1) << " a " << std::get<2>(tr2_mc1) << std::endl;
    } else {
        std::cout << "Mcode1 no se encuentra dentro de Transmision2" << std::endl;
    }

    std::cout << "\nParte 2: Busqueda de codigo espejeado" << std::endl;
    std::pair<int, int> result = buscador.manacher(transmission1);

    if (result.first == -1 && result.second == -1) {
        std::cout << "No se ha encontrado codigo espejeado" << std::endl;
    } else {
        printf("Se encontro codigo espejeado desde el indice %d hasta %d", result.first, result.second);
        std::cout << "\n";
    }
    
    std::pair<int, int> result2 = buscador.manacher(transmission2);

    if (result2.first == -1 && result2.second == -1) {
        std::cout << "No se ha encontrado codigo espejeado" << std::endl;
    } else {
        printf("Se encontro codigo espejeado desde el indice %d hasta %d", result2.first, result2.second);
        std::cout << "\n";
    }

    std::cout << "\n";

    std::cout << "\nParte 3: Coincidencia en transmiciones" << std::endl;

    buscador.lcs(transmission1, transmission2);

    return 0;
}
