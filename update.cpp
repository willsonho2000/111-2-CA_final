#include <iostream>
#include <math.h>
#include "octree.h"
#include "treewalk.h"

void tree_update( Octree* tree, double timestep, double** g ) {
    
    for ( int i = 0; i < tree->NumNodes; i++ ) {
        
        // walk through every particle in the tree
        Octree   *nod_target = tree->partree_arr[i];
        Particle *par_target = nod_target->par;

        // update the positions and velocity of the particle
        for ( int j = 0; j < 3; j++ ) {
            par_target->pos[j] += par_target->vel[i]*timestep + 0.5*g[i][j]*pow( timestep, 2 );
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
        
        std::cout << par_pos[0] << " " << par_pos[1] << " " << par_pos[2] << " " << reinsert << "\n";

        if ( reinsert ) {
            nod_target->par = nullptr;
            
            int i_octant = tree->FindQuad( par_target->pos, tree->Coordinates );
            tree->Insert( par_target, i_octant );
        }
    }
}