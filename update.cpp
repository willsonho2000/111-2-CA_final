#include <iostream>
#include <math.h>
#include "octree.h"
#include "treewalk.h"

void remove_node (Octree* pre_node, int to_octant, Octree* tree, double* node_pos) {
    
    // Don't proceed for nullptr
    if (tree == nullptr) return;

    // Find where the node is
    int octant = tree->FindQuad(node_pos, tree->Coordinates);
    Octree *target = tree->children[octant];
    
    int remove = 0;
    for (int i = 0; i < 3; i++) {
        if (target->Coordinates[i] == node_pos[i])
            remove = 1;
    }

    // if find the node, free the memory
    if (remove == 1) {
        delete[] target->Coordinates;
        target->children.clear();

        delete target;
        tree->children[octant] = nullptr;
    }
    else {  // or, keep going to the next level
        remove_node(tree, octant, tree->children[octant], node_pos);
    }
    
    // remove unecessary node
    int count = 0;
    for (int i = 0; i < 8; i++) {
        if (tree->children[i] != nullptr) count++;
    }

    if (count == 0) {

        // remove
        delete[] tree->Coordinates;
        for( int i = 0; i < 3; i++ ){
            delete[] tree->Quadrupoles[i];
        }
        delete[] tree->Quadrupoles;
        delete[] tree->com;
        tree->children.clear();

        delete tree;
        pre_node->children[to_octant] = nullptr;
    }

    return;
}

void tree_update (Octree* tree, double timestep, double** g) {
    
    for ( int i = 0; i < tree->NumNodes; i++ ) {
        
        // walk through every particle in the tree
        Octree   *nod_target = tree->partree_arr[i];
        Particle *par_target = nod_target->par;

        // update the positions and velocity of the particle
        for ( int j = 0; j < 3; j++ ) {
            par_target->pos[j] += par_target->vel[j]*timestep + 0.5*g[i][j]*pow( timestep, 2 );
            par_target->vel[j] += g[i][j]*timestep;
        }
        
        // check whether particle run outside the grid
        int reinsert = 0;
        double *par_pos  = par_target->pos;
        double *tree_pos = nod_target->Coordinates;
        double  size     = nod_target->Sizes;


        for ( int j = 0; j < 3; j++ ) {
            if ( par_pos[j] > tree_pos[j] + size || par_pos[j] < tree_pos[j] - size ) 
                reinsert = 1;
        }
        
        if ( reinsert ) {
            nod_target->par = nullptr;
            
            int i_octant = tree->FindQuad( par_target->pos, tree->Coordinates );
            tree->Insert( par_target, i_octant );

            remove_node(nullptr, 0, tree, nod_target->Coordinates);
        }
    }

    // Calculate the quadrupole of the new tree
    double mass, com[3], hmax; 
    ComputeMoments( &mass, com, &hmax, tree );
}
