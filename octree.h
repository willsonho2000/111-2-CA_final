#include <iostream>
#include <vector>

using namespace std;

struct Node {
    double pos[3];
    double mass;
    double softening;

    Node();
    Node(double* position, double m, double soft);
    Node(double a, double b, double c, double m, double soft);
    Node(int* position, double m, double soft);
    Node(int a, int b, int c, double m, double soft);
};

class Octree {
public:
    Node node;
    vector<Octree*> children;
    
    int NumParticles;
    int NumNodes;
    int Sizes;

    double* Coordinate;
    int* Deltas;
    int* Masses;

    Octree();   // initialization
    Octree(Node* root_node, double* coord, double size);

    void Insert(Node* new_node, int octant);
    // maybe can just call it once
    void BuildTree(double** points, double* masses, double* softenings, bool morton_order, bool quadrupole);
    // decide wich quad the node will be inserted to
    int FindQuad(double* pos, double* ref);
};