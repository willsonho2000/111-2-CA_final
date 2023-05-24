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

Node::Node( double* position, double m, double soft ) {
    for ( int i = 0; i < 3; i++ ) pos[i] = position[i];
    mass = m;
    softening = soft;
}

Node::Node( double a, double b, double c, double m, double soft ) {
    pos[0] = a;
    pos[1] = b;
    pos[2] = c;
    mass = m;
    softening = soft;
}

Node::Node( int* position, double m, double soft ) {
    for ( int i = 0; i < 3; i++ ) pos[i] = position[i];
    mass = m;
    softening = soft;
}

Octree::Octree() {
    this->node = nullptr;
    // Assigning null to the children
    children.assign( 8, nullptr );

    this->NumNodes = 0;
    this->Sizes = 0.;
    this->Coordinate = new double[3];
    this->root = this;
}

Octree::Octree( Node* root_node, Octree* root_ptr ) {
    this->node = root_node;
    // Assigning null to the children
    children.assign( 8, nullptr );

    this->root = root_ptr;
}

void Octree::Insert( Node* new_node, int octant ) {

    // check if there is a pre-existing particle among node's children
    if ( this->children[octant] != nullptr ) {
        // it's a particle
        if ( this->children[octant]->node != nullptr ) {
            
            // record the pre-existing particle's pointer
            Node* child_node = this->children[octant]->node;
            
            // EXCEPTION: if the pre-existing node is at the same coordinate, purturb the particle's position
            bool same_cood = true;
            for ( int i = 0; i < 3; i++ ) 
                if ( new_node->pos[i] != child_node->pos[i] ) same_cood = false;
            
            // restart the tree traversal (back to the root of the tree)
            if ( same_cood ) {
                srand (time(NULL));
                for ( int i = 0; i < 3; i++ ) {
                    new_node->pos[i] *= exp( 3e-16 * ((double) rand() / (RAND_MAX) - 0.5) );
                }

                // insert from the root 
                int new_octant = FindQuad( new_node->pos, this->root->Coordinate );
                this->root->Insert( new_node, new_octant );
            } // end exception

            // turn the node to a tree
            this->children[octant]->node = nullptr;

            // set the tree coordinate and sizes
            this->children[octant]->Sizes = this->Sizes/2.;
            this->children[octant]->Coordinate = new double[3];

            // set the value of coordinate
            double* cur_coord = this->Coordinate;
            double* child_coor = this->children[octant]->Coordinate;
            for ( int i = 0; i < 3; i++ ) 
                child_coor[i] = cur_coord[i] + this->Sizes*0.25*(double)octant_offset[octant][i];

            // put the pre-existing particle into the new tree
            int child_octant = FindQuad( child_node->pos, child_coor );
            this->children[octant]->children[child_octant] = new Octree(child_node, this->root);
            
            // insert new node in the new tree again
            int new_octant = FindQuad( new_node->pos, child_coor );
            this->children[octant]->Insert( new_node, new_octant );

        } // it's a tree, then pass the node to that tree
        else {
            int next_octant = FindQuad( new_node->pos, this->children[octant]->Coordinate );
            this->children[octant]->Insert( new_node, next_octant );
        }
    }
    else {

        // the child doesn't exist, so let the particle be the child
        this->children[octant] = new Octree( new_node, this->root );
    }

    this->NumNodes += 1;    // everytime we insert a node, parent's NumNodes + 1
}

void Octree::BuildTree( double** points, double* masses, double* softenings, bool morton_order=true, bool quadrapole=false ) {
    
    // record the max and the min of x, y, z
    int NumParticles = this->NumNodes = sizeof( points );
    
    for ( int i = 0; i < 3; i++ ) {

        // record the max and the min of the given axis
        double i_min, i_max;
        i_min = i_max = points[0][i];
        for ( int j = 0; j < NumParticles; j++ ) {
            double i_value = points[j][i];

            if ( i_max > i_value ) i_max = i_value;
            if ( i_min < i_value ) i_min = i_value;
        }

        // save properties of the grid size and the center of the position
        this->Coordinate[i] = 0.5 * ( i_max + i_min );
        if ( i_max - i_min > this->Sizes ) 
            this->Sizes = i_max - i_min;
    }

    // store the data into a new node and insert it into the tree
    for ( int i = 0; i < NumParticles; i++ ) {
        double* pos = points[i];

        // delcare a new node then insert it
        Node* i_node = new Node( pos, masses[i], softenings[i] );
        
        int i_octant = FindQuad( pos, this->Coordinate );
        this->Insert( i_node, i_octant );
    }
}

int Octree::FindQuad( double* pos, double* ref ) {
    
    int octant = 0;

    for ( int i = 0; i < 3; i++ )
        if ( pos[i] > ref[i] ) octant += 1 << i;
    
    return octant;
}
