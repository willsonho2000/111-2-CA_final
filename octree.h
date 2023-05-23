#include <iostream>
#include <vector>

using namespace std;

// double** octant_offset;

struct Node {
    double x;
    double y;
    double z;
    double mass;
    double softening;

    Node();
    Node(double a, double b, double c, double m, double soft);
};

class Octree {
    Node node;
    vector<Octree*> children;

public:
    vector<Octree*> BuildTree(double** points, double* masses, double* softenings);
};