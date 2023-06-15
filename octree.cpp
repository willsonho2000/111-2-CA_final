#include <iostream>
#include <vector>
#include "octree.h"
#include "treewalk.h"

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
    vel[0] = vel[1] = vel[2] = 0.;
    mass = -1.;
    softening = -1.;
}


Particle::Particle( double* position, double m, double soft ) {
    for ( int i = 0; i < 3; i++ ) pos[i] = position[i];
    vel[0] = vel[1] = vel[2] = 0.;
    mass = m;
    softening = soft;
}

Particle::Particle( double a, double b, double c, double m, double soft ) {
    pos[0] = a;
    pos[1] = b;
    pos[2] = c;
    vel[0] = vel[1] = vel[2] = 0.;
    mass = m;
    softening = soft;
}

Particle::Particle( int* position, double m, double soft ) {
    for ( int i = 0; i < 3; i++ ) pos[i] = position[i];
    vel[0] = vel[1] = vel[2] = 0.;
    mass = m;
    softening = soft;
}

// initialization
Octree::Octree( int N, double** points, double* masses, double* softening ) {

    this->par = nullptr;
    // Assign nullptr to the children
    children.assign( 8, nullptr );

    this->root = this;
    this->Coordinates = new double[3]{0};
    this->com         = new double[3]{0};  
    this->Quadrupoles = new double*[3];
    for ( int i = 0; i < 3; i++ ) this->Quadrupoles[i] = new double[3]{0.};
    this->NumNodes = N;
    this->Sizes = 0.;
    this->Softenings = 0.;
    this->Masses = 0.;
    // Assign nullptr to the particle array (only do this for root node)
    partree_arr.assign( N, nullptr );

    this->BuildTree( points, masses, softening );
    double mass, com[3], hmax; 
    ComputeMoments( &mass, com, &hmax, this );
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

                return;
            } // end exception

            // turn the node to a tree
            this->children[octant]->par = nullptr;

            // set the tree Coordinates, quadrupoles and sizes
            this->children[octant]->Sizes = this->Sizes/2.;
            this->children[octant]->Coordinates = new double[3]{0.};
            this->children[octant]->com = new double[3]{0.};
            this->children[octant]->Quadrupoles = new double*[3];
            for ( int i = 0; i < 3; i++ ) this->children[octant]->Quadrupoles[i] = new double[3]{0.};

            // set the value of Coordinates
            double* cur_coord = this->Coordinates;
            double* child_coor = this->children[octant]->Coordinates;
            for ( int i = 0; i < 3; i++ ) {
                child_coor[i] = cur_coord[i] + this->Sizes*0.25*( double )octant_offset[octant][i];
            }

            // put the pre-existing particle into the new tree
            int child_octant = FindQuad( child_par->pos, child_coor );
            this->children[octant]->children[child_octant] = new Octree( child_par, this->root );
            this->root->partree_arr[child_par->index] = this->children[octant]->children[child_octant];
            
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
        this->root->partree_arr[new_par->index] = this->children[octant];
    }
}

void Octree::BuildTree( double** points, double* masses, double* softenings ) {
    
    // record the max and the min of x, y, z
    for ( int i = 0; i < 3; i++ ) {

        // record the max and the min of the given axis
        double i_min, i_max;
        i_min = i_max = points[0][i];
        for ( int j = 0; j < this->NumNodes; j++ ) {
            double i_value = points[j][i];

            i_max = std::max( i_max, i_value );
            i_min = std::min( i_min, i_value );
        }

        // save properties of the grid size and the center of the position
        this->Coordinates[i] = 0.5 * ( i_max + i_min );
        if ( i_max - i_min > this->Sizes ) 
            this->Sizes = i_max - i_min;
    }
    // store the data into a new node and insert it into the tree
    for ( int i = 0; i < this->NumNodes; i++ ) {
        double* pos = points[i];

        // delcare a new node then insert it
        Particle* i_par = new Particle( pos, masses[i], softenings[i] );
        i_par->index = i;
        // this->par_arr[i] = i_par;

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