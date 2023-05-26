#include <iostream>
// #include "octree.h"

using namespace std;

int sum(int** arr) {
    int sum = 0;
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            sum += arr[i][j];
        }
    }

    return sum;
}

int main() {
    // int** a;
    // a = new int* [5];

    // for (int i = 0; i < 5; i++) {
    //     a[i] = new int [5];
    //     for (int j = 0; j < 5; j++) {
    //         a[i][j] = i+j;
    // }

    // cout << sum(a) << "\n";

    // int set[] = { 1,2,3 };
    // Particle a;
    // a = Particle( set, -1, -1. );
    // cout << a.pos[0] << "\n";
    // Octree *b = new Octree();

    // cout << *b->Coordinates << "\n";
    // cout << (b->node == nullptr) << "\n";

    std::cout << 2.4*2 << "\n";
    srand (time(NULL));
    // int random = rand();
    for (int i = 0; i < 5; i++)
        std::cout << ((double) rand() / (RAND_MAX) - 0.5) << " ";

    std::cout << "\n";

    double* b;
    b = new double[3]{0};

    std::cout << b[0] << "\n";

    return 0;
}