#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

class PrefixTree {
protected:
    struct Node {
        std::unordered_map<char, Node*> children;
        bool isEnd;
        char val;

        Node(bool isEnd, char val) : isEnd(isEnd), val(val) {}
    };

    Node* root;

    // Helper method for KMP Algorithm to compute the longest prefix suffix
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
    PrefixTree() : root(new Node(false, ' ')) {}

    // Insert a word into the trie
    bool insert(std::string word) {
        Node* current = root;
        for (char ch : word) {
            if (!current->children.count(ch)) {
                current->children[ch] = new Node(false, ch);
            }
            current = current->children[ch];
        }
        current->isEnd = true;
        return true;
    }

    // KMP Search
    bool KMP_Search(const std::string& text, const std::string& pattern) {
        std::vector<int> kmpTable = computeKMPTable(pattern);
        int i = 0; // index for text[]
        int j = 0; // index for pattern[]

        while (i < text.size()) {
            if (pattern[j] == text[i]) {
                j++;
                i++;
            }

            if (j == pattern.size()) {
                return true; // found pattern
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

    // Method to read file content into a string
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
    PrefixTree trie;
    std::string text = trie.readFileIntoString("transmission1.txt");
    trie.insert(text);

    std::string pattern = trie.readFileIntoString("mcode1.txt");
    bool found = trie.KMP_Search(text, pattern);
    std::cout << (found ? "Pattern found" : "Pattern not found") << std::endl;

    return 0;
}
