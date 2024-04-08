#include <iostream>
#include <unordered_map>
#include <vector>

class PrefixTree {
protected:
    struct Node {
        std::unordered_map<char, Node*> children;
        bool isEnd;
        char val;

        Node(bool isEnd, char val);
    };

    void insertHelper(std::vector<char> wordAsVector, Node* current);

    bool findHelper(std::vector<char> wordAsVector, Node* current);
public:
    Node* root;
    // Set root to nullptr
    PrefixTree();

    // Insert a word into the trie
    bool insert(std::string word);

    // Look for a word in the trie
    bool find(std::string word);
};