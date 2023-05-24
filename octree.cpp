#include <iostream>
#include <vector>
#include "octree.h"

using namespace std;

int octant_offset[8][3] =  {{-1,-1,-1},
                            {1,-1,-1},
                            {-1,1,-1},
                            {1,1,-1},
                            {-1,-1,1},
                            {1,-1,1},
                            {-1,1,1},
                            {1,1,1}};

Node::Node() {
    pos[0] = pos[1] = pos[2] = 0.;
    mass = -1.;
    softening = -1.;
}

Node::Node(double* position, double m, double soft) {
    for (int i = 0; i < 3; i++) pos[i] = position[i];
    mass = m;
    softening = soft;
}

Node::Node(double a, double b, double c, double m, double soft) {
    pos[0] = a;
    pos[1] = b;
    pos[2] = c;
    mass = m;
    softening = soft;
}

Node::Node(int* position, double m, double soft) {
    for (int i = 0; i < 3; i++) pos[i] = position[i];
    mass = m;
    softening = soft;
}

Node::Node(int a, int b, int c, double m, double soft) {
    pos[0] = a;
    pos[1] = b;
    pos[2] = c;
    mass = m;
    softening = soft;
}

// initialization
Octree::Octree() {
    node = Node();

    this->Sizes = 0.;
    this->Coordinate = new double[3];
    // this->Deltas = new int[8];
    // this->Masses = new int[8];
}

Octree::Octree(Node* root_node, double* coord, double size) {
    this->node = *root_node;

    this->Sizes = size;
    this->Coordinate = coord;
    // this->Deltas = new int[8];
    // this->Masses = new int[8];
}

void Octree::Insert(Node* new_node, int octant) {
    
    if (this->children[octant] != nullptr) {
        int next_octant = FindQuad(new_node->pos, this->children[octant]->Coordinate);
        this->children[octant]->Insert(new_node, next_octant);
    }
    else {
        // determine the center of the new node
        double* child_coordinate = new double[3];
        double* cur_coord = this->Coordinate;
        
        for (int i = 0; i < 3; i++) 
            child_coordinate[i] = cur_coord[i] + (double)Sizes*0.25*(double)octant_offset[octant][i];

        this->children[octant] = new Octree(new_node, child_coordinate, this->Sizes/2.);
    }
}

void Octree::BuildTree(double** points, double* masses, double* softenings, bool morton_order=true, bool quadrapole=false) {
    
    // record the max and the min of x, y, z
    this->NumParticles = sizeof(points);
    
    for (int i = 0; i < 3; i++) {
        // record the max and the min of the given axis
        double i_min, i_max;
        i_min = i_max = points[0][i];
        for (int j = 0; j < this->NumParticles; j++) {
            double i_value = points[j][i];

            if (i_max > i_value) i_max = i_value;
            if (i_min < i_value) i_min = i_value;
        }

        // save properties of the grid size and the center of the position
        this->Coordinate[i] = 0.5 * (i_max + i_min);
        if (i_max - i_min > this->Sizes) 
            this->Sizes = i_max - i_min;
    }

    // store the data into a new node and insert it into the tree
    for (int i = 0; i < this->NumParticles; i++) {
        double* pos = points[i];

        // delcare a new node then insert it
        Node* i_node = new Node(pos, masses[i], softenings[i]);
        
        int i_octant = FindQuad(pos, this->Coordinate);
        this->Insert(i_node, i_octant);
    }
}

int Octree::FindQuad(double* pos, double* ref) {
    int octant = 0;

    for (int i = 0; i < 3; i++)
        if (pos[i] > ref[i]) octant += 1 << i;
    
    return octant;
}
