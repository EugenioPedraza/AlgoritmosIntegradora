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
    if (wordAsVector.empty()) {
        current->isEnd = true;
        return;
    }
    else {
        auto it = current->children.find(wordAsVector[0]);

        if (it != current->children.end()){
            // Element found
            std::vector<char>newVect(wordAsVector.begin() + 1, wordAsVector.end());

            return insertHelper(newVect, current->children[wordAsVector[0]]);
        }
        else {
            // Element was not found 
            // Create a new node and set it up
            Node* newNode = new Node(false, wordAsVector[0]);

            // Append the new node into the current node's children
            current->children[wordAsVector[0]] = newNode;

            // Recursively call the function with the new node and without the newly added character
            std::vector<char>newVect(wordAsVector.begin() + 1, wordAsVector.end());
            return insertHelper(newVect, newNode);

        }
    }
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
    if (wordAsVector.empty()) {
        // Base case, the word wasn't found
        std::cout << "Character -> " << current->val << " found" << std::endl;
        return current->isEnd;
    }
    else {
        auto it = current->children.find(wordAsVector[0]);

        if (it != current->children.end()){
            // Element found
            std::vector<char>newVect(wordAsVector.begin() + 1, wordAsVector.end());
            if (current->val != ' '){
                std::cout << "Character -> " << current->val << " found" << std::endl;
            }
            return findHelper(newVect, current->children[wordAsVector[0]]);
        }

        return false;
    }

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

    std::cout << "Casos de prueba \n Buscando palabras...\n";

    std::cout << (trie.find(mcode1) ? "Transmission1 Encontrado\n" : "No encontrado\n");

    bool isActive = true;
    PrefixTree userTrie;


    return 0;
}