#include <iostream>
#include <vector>

struct Particle {
    double pos[3];
    double mass;
    double softening;

    Particle();
    Particle( double* position, double m, double soft );
    Particle( double a, double b, double c, double m, double soft );
    Particle( int* position, double m, double soft );
};

class Octree {
public:
    Particle* par;
    std::vector<Octree*> children;
    
    int NumNodes;           // how many particles the grid contain
    double Sizes;           // the size of the grid
    double Deltas;          // set deltas
    double Softenings;      // set softentings
    double Masses;          // set masses

    Octree* root;           // the root of the tree
    double* Coordinates;    // the center of the grid
    double** Quadrupoles;   // set quadrapoles

    // initialization an empty tree
    Octree( int N, double** points, double* masses, double* softening );
    Octree( Particle* root_par, Octree* root_ptr );    // initialize a new particle

    void Insert( Particle* new_par, int octant );
    void BuildTree( double** points, double* masses, double* softenings ); // maybe can just call it once
    int FindQuad( double* pos, double* ref );  // decide which quad the particle will be inserted to
};
