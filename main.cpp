#include "trie.h"
#include <fstream>
#include <iostream>
#include <string>
#include <iostream>
#include <sstream>

// Constructors
PrefixTree::Node::Node(bool isEnd, char val) {
    this->isEnd = isEnd;
    this->val = val;
}

PrefixTree::PrefixTree() {
    root = new Node(false, ' ');
}

// Insert words into trie
bool PrefixTree::insert(std::string word) {
    if (word.empty()) {
        return false;
    }

    std::vector<char> wordVect(word.begin(), word.end());
    insertHelper(wordVect, root);

    return true;
}

// Recursive insert helper
void PrefixTree::insertHelper(std::vector<char> wordAsVector, Node* current) {
    for (int i = 0; i <  wordAsVector.size(); i++) {
        auto it = current->children.find(wordAsVector[i]);

        if (it != current->children.end()){
            // Element found
            current = it->second; // Move to the next node in the tree
        }
        else {
            // Element was not found 
            // Create a new node and set it up
            Node* newNode = new Node(false, wordAsVector[i]);

            // Append the new node into the current node's children
            current->children[wordAsVector[i]] = newNode;

            // Move to the next node in the tree
            current = newNode;
        }
    }

    current->isEnd = true; // Mark the last node as the end of a word
}

// Find words
bool PrefixTree::find(std::string word) {
    if (word.empty()) {
        return false;
    }

    std::vector<char> wordVect(word.begin(), word.end());
    return findHelper(wordVect, root);
} 

// Recursive find helper 
bool PrefixTree::findHelper(std::vector<char> wordAsVector, Node* current) {
    while (!wordAsVector.empty()) {
        auto it = current->children.find(wordAsVector[0]);

        if (it != current->children.end()){
            // Element found
            std::vector<char> newVect(wordAsVector.begin() + 1, wordAsVector.end());
    
            current = current->children[wordAsVector[0]];
            wordAsVector = newVect;
        }
        else {
            return false;
        }
    }

    // Base case, the word wasn't found
    return current->isEnd;
}

std::string readFileIntoString(const std::string& path) {
    std::ifstream inputFile(path);
    if (!inputFile.is_open()){
        std::cerr << "Error abriendo archivo " << path <<std::endl;
        return "";
    }
    std::stringstream buffer;
    buffer << inputFile.rdbuf();
    inputFile.close();
    return buffer.str();
}

int main() {
    PrefixTree trie;
    std::string mcode1, transmission1, transmission2;

    mcode1 = readFileIntoString("mcode1.txt");
    transmission1 = readFileIntoString("transmission1.txt");
    transmission2 = readFileIntoString("transmission2.txt"); 

    trie.insert(transmission1);
    trie.insert(transmission2);

    std::cout << trie.find(transmission2);

    bool isActive = true;
    PrefixTree userTrie;


    return 0;
}