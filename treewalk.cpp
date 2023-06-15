#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <omp.h>
#include "octree.h"
#include "treewalk.h"
using namespace std;

void ComputeMoments( double* mass, double com[3], double* hmax, Octree* tree )
{
    // #####################
    // Note:
    // 1. quad is intitially 0
    // #####################
    if ( tree == nullptr ) return;

    bool IsParticle = true; // checking whether the current node is a particle
    if ( tree->par == nullptr )   IsParticle = false;

    double** quad;
    quad = new double*[3]; // empty
    for (int i = 0; i < 3; i++)  quad[i] = new double[3]{0.};

    if ( IsParticle )
    {
        *mass = tree->par->mass;
        for ( int i = 0; i < 3; i++ ) com[i] = tree->par->pos[i];
        *hmax = tree->par->softening;

        return;
    }
    else // it is a node
    {
        double hmax0   = 0.;                // some properties for the node
        double m0      = 0.;                // total mass
        double com0[3]{0.};                 // used to calculate COM of the node
        double comi[3];                     // some properties for the child

        for ( int octant=0; octant<8; octant++ ) // open the node to calculate the total mass and COM position
        {
            // some properties for the child
            double hmaxi;
            double massi;

            if ( tree->children[octant] == nullptr ) continue;
            ComputeMoments( &massi, comi, &hmaxi, tree->children[octant] );

            hmax0 = max( hmax0, hmaxi );
            m0 += massi;
            for ( int i = 0; i < 3; i++ ) com0[i] += massi*comi[i];
        } // for (int octant=0; octant<8; octant++)

        // compute the COM
        for ( int i = 0; i < 3; i++ ) com0[i] = com0[i]/m0;

        for ( int octant=0; octant<8; octant++ ) // open the node to calculate quadrapoles from children
        {
            if ( tree->children[octant] == nullptr ) continue;
            
            double ri[3]{0.};
            double r2 = (double)0;

            double*  comi;
            double** quadi;
            double mi;

            if ( tree->children[octant]->par != nullptr ) {
                quadi = quad;
                comi  = tree->children[octant]->par->pos;
                mi    = tree->children[octant]->par->mass;
            }
            else {
                quadi = tree->children[octant]->Quadrupoles;
                comi  = tree->children[octant]->com;
                mi    = tree->children[octant]->Masses;
            }

            for (int i = 0; i < 3; i++ ) ri[i] = comi[i] - com0[i];
            for (int i = 0; i < 3; i++ ) r2 += ri[i]*ri[i];

            for (int k = 0; k < 3; k++ )
            for (int l = 0; l < 3; l++ )
            {
                if ( tree->children[octant]->par != nullptr )
                {
                    quad[k][l] += mi*3*ri[k]*ri[l];
                    if ( k==l ) quad[k][l] -= mi*r2;
                }
                else // if the child is a node, just propagate quadi up to the parent
                {
                    quad[k][l] += quadi[k][l];
                }
            } // l, k
        } // for (int octant=0; octant<8; octant++)

        double delta = (double)0;
        for ( int dim = 0; dim < 3; dim++ )
        {
            double dx = com0[dim] - tree->com[dim];
            delta += pow( dx, 2 );
        } // for (int dim=0; dim<3; dim++)

        // update tree properties
        tree->Masses      = m0;
        for ( int i = 0; i < 3; i++ ) tree->com[i] = com0[i];
        tree->Softenings  = hmax0;
        double** quadi = tree->Quadrupoles;
        tree->Quadrupoles = quad;
        tree->Deltas      = sqrt( delta );


        // return the properties to parent node
        *mass       = tree->Masses;
        for ( int i = 0; i < 3; i++ ) com[i] = com0[i];
        *hmax       = tree->Softenings;

        if ( tree->par == nullptr ) {
        
            // free memory
            for( int i = 0; i < 3; i++ ){
                delete[] quadi[i];
            }
        }

        delete[] quadi;

        return;
    }
}

double PotentialKernel(double r, double h)
{
    if ( h==0.0 )
    {
        return -1/r;
    }

    double hinv = 1/h;
    double q = r*hinv;

    if ( q<=0.5 )
    {
        return (-2.8 + q*q*(5.33333333333333333 + q*q*(6.4*q - 9.6))) * hinv;
    }
    else if ( q <= 1 )
    {
        return (-3.2 + 0.066666666666666666666 / q + q*q*(10.666666666666666666666 +  q*(-16.0 + q*(9.6 - 2.1333333333333333333333 * q)))) * hinv;
    }
    else
    {
        return -1./r;
    }
}

double ForceKernel(double r, double h)
{
    if ( r>h )
    return 1./(r*r*r);

    double hinv = 1./h;
    double q = r*hinv;

    if ( q<=0.5 )
    {
        return (10.666666666666666666 + q*q*(-38.4 + 32.*q))*hinv*hinv*hinv;
    }
    else
    {
        return (21.333333333333 - 48.0 * q + 38.4 * q * q - 10.666666666667 * q * q * q - 0.066666666667 / (q * q * q))*hinv*hinv*hinv;
    }
}

double PotentialWalk_quad(double* pos, Octree* tree, double theta, double softening)
{
    if ( tree == nullptr ) return 0.0;

    double phi = 0;
    double r5inv;

    bool IsParticle = true; // checking whether the current node is a particle
    if ( tree->par == nullptr )   IsParticle = false;

    double r=0, dx[3] = {};
    for (int k=0; k<3; k++)
    {   
        if ( IsParticle )
        {
            dx[k] = tree->par->pos[k] - pos[k];
        }
        else
        {
            dx[k] = tree->com[k] - pos[k];
        }
        r+=dx[k]*dx[k];
    }
    r = sqrt(r); // distance between observer and the node

    double h = fmax(tree->Softenings, softening); // softening length

    if ( IsParticle ) // it is a particle
    {
        if ( r > 0 ) // itself, we don't calculate self-gravity
        {
            if ( r < h)
            {
                phi += tree->par->mass * PotentialKernel(r,h); 
            }
            else
            {
                phi -= tree->par->mass / r;
            } // if ( r < h)
        } // if ( r > 0 )
    } // if ( IsParticle )
    else if ( r > fmax(tree->Sizes/theta + tree->Deltas, h+tree->Sizes*0.6+tree->Deltas) ) // satisfy opening criteria
    {
        phi -= tree->Masses/r;
        r5inv = 1 / pow(r, 5);

        for (int k=0; k<3; k++)
        for (int l=0; l<3; l++)
        phi -= 0.5 * dx[k] * tree->Quadrupoles[k][l] * dx[l] * r5inv;
    } // else if ( r > fmax(tree->Sizes/theta + tree->Deltas, h+tree->Sizes*0.6+tree->Deltas) )
    else
    {
        for (int octant=0; octant<8; octant++) // open the node
        {
            phi += PotentialWalk_quad(pos, tree->children[octant], theta, softening);
        }
    } // else

    return phi;
}

double* AccelWalk_quad(double* pos, Octree* tree, double theta, double softening)
{
    double* g = new double[3]{0.};
    if ( tree == nullptr ) return g;

    double r, r5inv, fac;

    bool IsParticle = true; // checking whether the current node is a particle
    if ( tree->par == nullptr )   IsParticle = false;

    double r2=0, dx[3] = {};
    for (int k=0; k<3; k++)
    {   
        if ( IsParticle )
        {
            dx[k] = tree->par->pos[k] - pos[k];
        }
        else
        {
            dx[k] = tree->com[k] - pos[k];
        }
        r2+=dx[k]*dx[k];
    }
    r = sqrt(r2); // distance between observer and the node

    double h = fmax(tree->Softenings, softening); // softening length

    if ( IsParticle ) // it is a particle
    {
        if ( r > 0 ) // itself, we don't calculate self-gravity
        {
            if ( r < h)
            {
                fac = tree->par->mass * ForceKernel(r,h); 
            }
            else
            {
                fac = tree->par->mass/(r*r2);
            } // if ( r < h)
        } // if ( r > 0 )
        for (int k=0; k<3; k++) g[k] += fac * dx[k]; // monopole
    } // if ( IsParticle )
    else if ( r > fmax(tree->Sizes/theta + tree->Deltas, h+tree->Sizes*0.6+tree->Deltas) ) // satisfy opening criteria
    {
        double quad_fac = 0.0;
        fac = tree->Masses/(r*r2);
        r5inv = 1 / (r2*r2*r);

        for (int k=0; k<3; k++)
        {
            g[k] += fac * dx[k]; // monopole
            for (int l=0; l<3; l++)
            quad_fac += tree->Quadrupoles[k][l] * dx[k] * dx[l];
        }
        quad_fac *= r5inv / r2;

        for (int k=0; k<3; k++)
        {
            g[k] += 2.5 * quad_fac * dx[k];
            for (int l=0; l<3; l++)
            g[k] -= tree->Quadrupoles[k][l] * dx[l] * r5inv;
        }
    } // else if ( r > fmax(tree->Sizes/theta + tree->Deltas, h+tree->Sizes*0.6+tree->Deltas) )
    else
    {
        for (int octant=0; octant<8; octant++) // open the node
        {
            double* g_node;
            g_node = AccelWalk_quad(pos, tree->children[octant], theta, softening);

            for (int k=0; k<3; k++) g[k] += g_node[k];
        }
    } // else

    return g;
}

double* PotentialTarget_tree(int Npar, double** pos_target, double* softening_target, Octree* tree, double G, double theta)
{
    double* result = new double[Npar];

    printf( "Number of threads = %d for potential calculation.\n", omp_get_max_threads() );

#   pragma omp parallel for
    for (int i=0; i<Npar; i++)
    result[i] = G*PotentialWalk_quad(pos_target[i], tree, theta, softening_target[i]);

    return result;
}

double** AccelTarget_tree(int Npar, double** pos_target, double* softening_target, Octree* tree, double G, double theta)
{
    double** result = new double*[Npar];
    for (int i = 0; i < Npar; i++)  result[i] = new double[3];

    // printf( "Number of threads = %d for acceleration calculation.\n", omp_get_max_threads() );

#   pragma omp parallel for
    for (int i=0; i<Npar; i++)
    {
        double* g_tree;
        g_tree = AccelWalk_quad(pos_target[i], tree, theta, softening_target[i]);

        for (int k=0; k<3; k++) result[i][k] = G*g_tree[k];
    }

    return result;
}
