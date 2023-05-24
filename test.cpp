#include <iostream>
#include "octree.h"

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

// struct Node {
//     double pos[3];
//     double mass;
//     double softening;

//     Node();
//     // Node(double a, double b, double c, double m, double soft);
// };

// Node::Node() {
//     pos[0] = pos[1] = pos[2] = -1.;
//     mass = -1.;
//     softening = -1.;
// }

int main() {
    // int** a;
    // a = new int* [5];

    // for (int i = 0; i < 5; i++) {
    //     a[i] = new int [5];
    //     for (int j = 0; j < 5; j++) {
    //         a[i][j] = i+j;
    //     }
    // }m

    // cout << sum(a) << "\n";

    int set[] = {1,2,3};
    Node a;
    a = Node(set, -1, -1.);
    cout << a.pos[0] << "\n";
    Octree *b = new Octree();

    cout << *b->Coordinate << "\n";
    cout << (b->node == nullptr) << "\n";

    cout << 2.4*2 << "\n";
    srand (time(NULL));
    int random = rand();
    cout << random%100 / (double) 100 << "\n";

    return 0;
}