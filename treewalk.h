#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <algorithm>
using namespace std;


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

double PotentialWalk_quad(double* pos, Octree* tree, double theta, double softening)
{
    double phi = 0;
    double quad[3][3], r5inv;

    double r=0, dx[3] = {};
    for (int k=0; k<3; k++)
    {
        dx[k] = tree->Coordinates[k] - pos[k];
        r+=dx[k]*dx[k];
    }
    r = sqrt(r); // distance between observer and the node

    double h = fmax(tree->Softenings, softening); // softening length

    bool IsParticle = true;
    for (int octant=0; octant<8; octant++) // check whether the current node is particle
    {
        if ( tree->children[octant] != nullptr )
        {
            IsParticle = false;
            break;
        }
    } // for (int octant=0; octant<8; octant++)

    if ( IsParticle ) // it is a particle
    {
        if ( r > 0 )
        {
            if ( r < h)
            {
                phi += tree->Masses * PotentialKernel(r,h); 
            }
            else
            {
                phi -= tree->Masses / r;
            } // if ( r < h)
        } // if ( r > 0 )
    } // if ( IsParticle )
    else if ( r > fmax(tree->Sizes/theta + tree->Deltas, h+tree->Sizes*0.6+tree->Deltas) ) // satisfy opening criteria
    {
        phi -= tree->Masses/r;
        copy(&tree->Quadrupoles[0][0], &tree->Quadrupoles[0][0] + 3 * 3, &quad[0][0]);
        r5inv = 1 / pow(r, 5);

        for (int k=0; k<3; k++)
        for (int l=0; l<3; l++)
        phi -= 0.5 * dx[k] * quad[k][l] * dx[l] * r5inv;
    } // else if ( r > fmax(tree->Sizes/theta + tree->Deltas, h+tree->Sizes*0.6+tree->Deltas) )
    else
    {
        for (int octant=0; octant<8; octant++) // open the node
        phi += PotentialWalk_quad(pos, tree->children[octant], theta, softening);
    } // else

    return phi;
}

double* PotentialTarget_tree(double** pos_target, double* softening_target, Octree* tree, int G, double theta)
{
    int N = sizeof(pos_target);
    double result[N];

    for (int i=0; i<N; i++)
    result[i] = G*PotentialWalk_quad(pos_target[i], tree, softening_target[i], theta);

    return result;
}