#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "res.h"


std::string Algorithm::manacher(std::string transmissionText){ 
    int N = transmissionText.size();
    if (N == 0)
        return "";
    N = 2*N + 1; // Position count
    std::vector<int> L(N, 0); // LPS Length Array
    L[0] = 0;
    L[1] = 1;
    int C = 1; // centerPosition
    int R = 2; // centerRightPosition
    int i = 0; // currentRightPosition
    int iMirror; // currentLeftPosition
    int maxLPSLength = 0;
    int maxLPSCenterPosition = 0;
    int start = -1;
    int end = -1;
    int diff = -1;

    for (i = 2; i < N; i++) {
        // get currentLeftPosition iMirror for currentRightPosition i
        iMirror = 2*C-i;
        L[i] = 0;
        diff = R - i;
        // If currentRightPosition i is within centerRightPosition R
        if (diff > 0)
            L[i] = std::min(L[iMirror], diff);

        // Attempt to expand palindrome centered at currentRightPosition i
        while (((i + L[i]) < N && (i - L[i]) > 0) && 
            (((i + L[i] + 1) % 2 == 0) || 
            (transmissionText[(i + L[i] + 1)/2] == transmissionText[(i - L[i] - 1)/2]))) {
            L[i]++;
        }

        if (L[i] > maxLPSLength) { // Track maxLPSLength
            maxLPSLength = L[i];
            maxLPSCenterPosition = i;
        }

        // If palindrome centered at currentRightPosition i 
        // expand beyond centerRightPosition R,
        // adjust centerPosition C based on expanded palindrome.
        if (i + L[i] > R) {
            C = i;
            R = i + L[i];
        }
    }
    start = (maxLPSCenterPosition - maxLPSLength)/2;
    end = start + maxLPSLength - 1;
    return transmissionText.substr(start, end-start+1);
}

class KMP_Searcher {
private:
    std::vector<int> computeKMPTable(const std::string& word) {
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

public:
    std::tuple<bool,int, int> KMP_Search(const std::string& text, const std::string& pattern) {
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

    
    std::string readFileIntoString(const std::string& path) {
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
};

int main() {
    KMP_Searcher buscador;
    Algorithm find;

    std::string transmission1 = buscador.readFileIntoString("transmission1.txt");
    std::string transmission2 = buscador.readFileIntoString("transmission2.txt");
    
    std::cout << find.manacher(transmission1);
}
