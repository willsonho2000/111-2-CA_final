#include <iostream>
#include <math.h>
#include "octree.h"
#include "treewalk.h"

void tree_update( Octree* tree, double timestep, double** g ) {
    
    for ( int i = 0; i < tree->NumNodes; i++ ) {
        
        // walk through every particle in the tree
        Particle* par_target = tree->par_arr[i];
        Octree*   nod_target = par_target->node;

        // update the positions and velocity of the particle
        for ( int j = 0; j < 3; j++ ) {
            par_target->pos[i] += par_target->vel[i]*timestep + 0.5*g[i][j]*pow( timestep, 2 );
            par_target->vel[i] += g[i][j]*timestep;
        }
        
        // check whether particle run outside the grid
        int reinsert = 0;
        if ( par_target )
    }
}

Octree* find_par( Octree* tree, int target ) {

}