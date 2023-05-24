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
    // }

    // cout << sum(a) << "\n";

    int set[] = {1,2,3};
    Node a;
    a = Node(-1, -2, -3, -1, -1);
    cout << a.pos[0] * 2.5 << "\n";
    Octree *b = new Octree();

    cout << sizeof(b->Coordinate) << "\n";
    cout << b->Coordinate[0] << "\n";
    cout << b->node.mass << "\n";

    return 0;
}