#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

class Algorithm {
    private:
    std::vector<int> computeKMPTable(const std::string& word);

    public:
    // Read file into string
    std::string readFileIntoString(const std::string& path);

    // KMP Search, part 1
    std::tuple<bool,int, int> KMP_Search(const std::string& text, const std::string& pattern);

    // Manacher's algorithm, part 2
    std::string manacher(std::string transmissionText);
};