#include <iostream>
#include <vector>

using namespace std;

struct Node {
    double pos[3];
    double mass;
    double softening;

    Node();
    Node( double* position, double m, double soft );
    Node( double a, double b, double c, double m, double soft );
    Node( int* position, double m, double soft );
};

class Octree {
public:
    Node* node;
    vector<Octree*> children;
    
    int NumNodes;       // how many particles the grid contain
    double Sizes;       // the size of the grid
    double* Coordinate; // the center of the grid
    Octree* root;       // the root of the tree

    Octree();                                   // initialization an empty tree
    Octree( Node* root_node, Octree* root_ptr );    // initialize a new particle

    void Insert( Node* new_node, int octant );
    void BuildTree( double** points, double* masses, double* softenings, bool morton_order, bool quadrupole ); // maybe can just call it once
    int FindQuad( double* pos, double* ref );  // decide which quad the node will be inserted to
};