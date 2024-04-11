#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

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
    bool KMP_Search(const std::string& text, const std::string& pattern) {
        std::vector<int> kmpTable = computeKMPTable(pattern);
        int i = 0;
        int j = 0; 

        while (i < text.size()) {
            if (pattern[j] == text[i]) {
                j++;
                i++;
            }

            if (j == pattern.size()) {
                return true;
                j = kmpTable[j - 1];
            } else if (i < text.size() && pattern[j] != text[i]) {
                if (j != 0) {
                    j = kmpTable[j - 1];
                } else {
                    i = i + 1;
                }
            }
        }

        return false; 
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
    std::string transmision1 = buscador.readFileIntoString("transmision1.txt");
    std::string transmision2 = buscador.readFileIntoString("transmision2.txt");
    std::string mcode1 = buscador.readFileIntoString("mcode1.txt");

    bool tr1_mc1 = buscador.KMP_Search(transmision1, mcode1);
    bool tr2_mc1 = buscador.KMP_Search(transmision2, mcode1);
    std::cout << (tr1_mc1 ? "Mcode1 se encuentra dentro de Transmision1" : "Mcode1 no se encuentra dentro de Transmision1") << std::endl;
    std::cout << (tr2_mc1 ? "Mcode1 se encuentra dentro de Transmision2" : "Mcode1 no se encuentra dentro de Transmision2") << std::endl;
    return 0;
}

