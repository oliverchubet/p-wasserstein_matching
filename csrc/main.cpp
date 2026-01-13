#include <iostream>
#include "data_loader.hpp" 

int main() {
    std::cout << "--- Reading data from files ---\n";
    const std::string inputFilenameA = "points_A.txt";
    const std::string inputFilenameB = "points_B.txt";
    Dataset A;
    Dataset B;
    if (readPointsFromFile(inputFilenameA, A)) {
        std::cout << "Loaded point set A.\n";
        //print_points(A);
    } else {
        std::cerr << "Data loading of point set A failed. Check file existence and permissions.\n";
        return 1;
    }
    if (readPointsFromFile(inputFilenameB, B)) {
        std::cout << "Loaded point set B.\n";
        //print_points(B);
    } else {
        std::cerr << "Data loading of point set B failed. Check file existence and permissions.\n";
        return 1;
    }


    return 0;
}
