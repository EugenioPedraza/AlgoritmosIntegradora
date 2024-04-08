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

    trie.insert(mcode1);

    // Insertar palabras

    // Buscar palabras
    std::cout << "Casos de prueba \n Buscando palabras...\n";

    std::cout << "hola: " << (trie.find(transmission1) ? "Transmission1 Encontrado\n" : "No encontrado\n");

    // Usuario 
    bool isActive = true;
    PrefixTree userTrie;

    while (isActive) {
        int action;
        std::string searchWord;

        std::cout << "Que desea hacer?" << std::endl;
        std::cout << "1) Agregar palabras" << std::endl;
        std::cout << "2) Buscar una palabra" << std::endl;
        std::cout << "3) Salir" << std::endl;

        std::cout << "Seleccion: ";
        std::cin >> action;

        if(std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }

        int quantity;
        int success = 0;
        int failure = 0;
        std::vector<std::string> words;

        switch (action)
        {
        case 1:
            std::cout << "Cuantas palabras desea insertar?: ";
            while(!(std::cin >> quantity)) {
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                std::cout << "Entrada invalida";
            }

            printf("Inserte %d palabras de una por una, una vez que haya escrito la palabra, presione enter para continuar, procure no utilizar espacios \n", quantity);

            for (int i = 0; i < quantity; i++) {
                printf("Palabras restantes: %d \n", quantity - i);
                std::string input;

                std::getline(std::cin >> std::ws, input);

                while(input.find(' ') != std::string::npos){
                    std::cout << "La entrada consiste de mas de 1 palabras o contiene un espacio, intente de nuevo" << std::endl;
                    std::cin >> input;
                } 

                words.push_back(input);      
            }

            for (int i = 0; i < words.size(); i++) {
                bool result = userTrie.insert(words[i]);

                if (result) {
                    success++;
                }
                else {
                    failure++;
                }
            }

            printf("%d palabras fueron agregadas con exito! \n", success);
            printf("%d palabras no fueron agregadas \n", failure);
            break;

        case 2:
            std::cout << "Ingrese la palabra que desea buscar: ";

            std::getline(std::cin >> std::ws, searchWord);

            while(searchWord.find(' ') != std::string::npos){
                std::cout << "La entrada consiste de mas de 1 palabras o contiene un espacio, intente de nuevo" << std::endl;
                std::cin >> searchWord;
            }

            if (userTrie.find(searchWord)) {
                std::cout << "La palabra '" << searchWord << "' fue encontrada en el trie!" << std::endl;
            } else {
                std::cout << "La palabra '" << searchWord << "' no fue encontrada en el trie." << std::endl;
            }

            break;

        case 3:
            isActive = false;
            break;

        default:
            printf("***La seleccion %d es invalida, intente de nuevo*** \n", action);
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            break;
        }
    }


    return 0;
}