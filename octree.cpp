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

Particle::Particle() {
    pos[0] = pos[1] = pos[2] = 0.;
    mass = -1.;
    softening = -1.;
}

Particle::Particle( double* position, double m, double soft ) {
    for ( int i = 0; i < 3; i++ ) pos[i] = position[i];
    mass = m;
    softening = soft;
}

Particle::Particle( double a, double b, double c, double m, double soft ) {
    pos[0] = a;
    pos[1] = b;
    pos[2] = c;
    mass = m;
    softening = soft;
}

Particle::Particle( int* position, double m, double soft ) {
    for ( int i = 0; i < 3; i++ ) pos[i] = position[i];
    mass = m;
    softening = soft;
}

Octree::Octree( double** points, double* masses, double* softening, bool morton_order=true, bool quadrupole=false ) {
    this->par = nullptr;
    // Assigning null to the children
    children.assign( 8, nullptr );

    this->root = this;
    this->Coordinates = new double[3]{0.};
    this->NumNodes = 0;
    this->Sizes = 0.;
    this->Softenings = 0.;
    this->Masses = 0.;

    this->BuildTree( points, masses, softening );
}

Octree::Octree( Particle* root_par, Octree* root_ptr ) {
    this->par = root_par;
    // Assigning null to the children
    children.assign( 8, nullptr );

    this->root = root_ptr;
}

void Octree::Insert( Particle* new_par, int octant ) {

    // check if there is a pre-existing particle among node's children
    if ( this->children[octant] != nullptr ) {
        // it's a particle
        if ( this->children[octant]->par != nullptr ) {
            
            // record the pre-existing particle's pointer
            Particle* child_par = this->children[octant]->par;
            
            // EXCEPTION: if the pre-existing node is at the same coordinate, purturb the particle's position
            bool same_cood = true;
            for ( int i = 0; i < 3; i++ ) 
                if ( new_par->pos[i] != child_par->pos[i] ) same_cood = false;
            
            // restart the tree traversal (back to the root of the tree)
            if ( same_cood ) {
                srand (time(NULL));
                for ( int i = 0; i < 3; i++ ) {
                    new_par->pos[i] *= exp( 3e-16 * ((double) rand() / (RAND_MAX) - 0.5) );
                }

                // insert from the root 
                int new_octant = FindQuad( new_par->pos, this->root->Coordinates );
                this->root->Insert( new_par, new_octant );
            } // end exception

            // turn the node to a tree
            this->children[octant]->par = nullptr;

            // set the tree Coordinates and sizes
            this->children[octant]->Sizes = this->Sizes/2.;
            this->children[octant]->Coordinates = new double[3];

            // set the value of Coordinates
            double* cur_coord = this->Coordinates;
            double* child_coor = this->children[octant]->Coordinates;
            for ( int i = 0; i < 3; i++ ) 
                child_coor[i] = cur_coord[i] + this->Sizes*0.25*( double )octant_offset[octant][i];

            // put the pre-existing particle into the new tree
            int child_octant = FindQuad( child_par->pos, child_coor );
            this->children[octant]->children[child_octant] = new Octree(child_par, this->root);
            
            // insert new node in the new tree again
            int new_octant = FindQuad( new_par->pos, child_coor );
            this->children[octant]->Insert( new_par, new_octant );

        } // it's a tree, then pass the node to that tree
        else {
            int next_octant = FindQuad( new_par->pos, this->children[octant]->Coordinates );
            this->children[octant]->Insert( new_par, next_octant );
        }
    }
    else {
        // the child doesn't exist, so let the particle be the child
        this->children[octant] = new Octree( new_par, this->root );
    }

    this->NumNodes += 1;    // everytime we insert a node, parent's NumNodes + 1
}

void Octree::BuildTree( double** points, double* masses, double* softenings ) {
    
    // initialization
    int NumParticles = this->NumNodes = sizeof( points );
    this->Quadrupoles = new double*[3];
    for ( int i = 0; i < NumParticles; i++ ) 
        Quadrupoles[i] = new double[3];
    
    // record the max and the min of x, y, z
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
        this->Coordinates[i] = 0.5 * ( i_max + i_min );
        if ( i_max - i_min > this->Sizes ) 
            this->Sizes = i_max - i_min;
    }

    // store the data into a new node and insert it into the tree
    for ( int i = 0; i < NumParticles; i++ ) {
        double* pos = points[i];

        // delcare a new node then insert it
        Particle* i_par = new Particle( pos, masses[i], softenings[i] );
        
        int i_octant = FindQuad( pos, this->Coordinates );
        this->Insert( i_par, i_octant );
    }
}

int Octree::FindQuad( double* pos, double* ref ) {
    
    int octant = 0;

    for ( int i = 0; i < 3; i++ )
        if ( pos[i] > ref[i] ) octant += 1 << i;
    
    return octant;
}