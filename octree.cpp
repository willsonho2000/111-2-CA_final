#include <iostream>
#include <vector>
#include "octree.h"

using namespace std;

// class Node {
// public:
//     int val;
//     Node* next;

//     Node(int x): val(x), next(nullptr) {}
//     Node(int x, Node* next): val(x), next(next) {}
// };

double** octant_offset = new double*[8];

Node::Node(): 
    x(-1), y(-1), z(-1), mass(-1), softening(-1)
    {}


Node::Node(double a, double b, double c, double m, double soft): 
    x(a), y(b), z(c), mass(m), softening(soft)
    {}

vector<Octree*> Octree::BuildTree(double** points, double* masses, double* softenings) {
    return children;
}
