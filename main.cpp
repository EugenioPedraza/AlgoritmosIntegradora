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
    std::string transmission1 = buscador.readFileIntoString("transmission1.txt");
    std::string transmission2 = buscador.readFileIntoString("transmission2.txt");
    std::string mcode1 = buscador.readFileIntoString("mcode1.txt");

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

    return 0;
}
