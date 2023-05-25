#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <algorithm>
using namespace std;

void ComputeMoments(double* mass, double* com, double** quad, double* hmax, Octree* tree)
{
    // #####################
    // Note:
    // 1. quad is intitially 0
    // #####################

    bool IsParticle = true; // checking whether the current node is a particle
    if ( tree->par == nullptr )   IsParticle = false;

    double* hmaxi; // some properties for the child
    double* massi;
    double*  comi  = new double[3]{0};
    double** quadi = new double*[3];
    for (int i = 0; i < 3; i++)  quadi[i] = new double[3]{0};

    quad = quadi; // empty

    if ( IsParticle )
    {
        mass = tree->par->mass;
        com  = tree->par->pos;
        hmax = tree->par->softening;
    }
    else // it is a node
    {
        double hmax0   = 0; // some properties for the node
        double m0      = 0; // total mass
        double com0[3] = {}; // used to calculate COM of the node

        for (int octant=0; octant<8; octant++) // open the node to calculate the total mass and COM position
        {
            ComputeMoments(massi, comi, quadi, hmaxi, tree->children[octant]);

            hmax0 = max(hmax0, *hmaxi);
            m0 += *massi;
            for (int i = 0; i < 3; i++) com0[i] += (*massi)*comi[i];
        } // for (int octant=0; octant<8; octant++)

        tree->Masses = mass = &m0;
        for (int i = 0; i < 3; i++) com[i] = com0[i]/m0;

        for (int octant=0; octant<8; octant++) // open the node to calculate quadrapoles from children
        {
            double* ri = new double[3];
            double  r2 = (double)0;

            comi  = tree->children[octant]->Coordinate;
            quadi = tree->children[octant]->Quadrapoles;

            for (int i = 0; i < 3; i++) ri[i] = comi[i] - com[i];
            for (int i = 0; i < 3; i++) r2 += ri[i]*ri[i];

            for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
            {
                quad[k][l] += quadi[k][l] + tree->children[octant]->Masses*3*ri[k]*ri[l];
                if ( k==l ) quad[k][l] -= tree->children[octant]->Masses*r2;
            } // l, k

            delete[] ri;
        } // for (int octant=0; octant<8; octant++)

        tree->Quadrapoles = quad;

        double delta = (double)0;
        for (int dim=0; dim<3; dim++)
        {
            double dx = com[dim] - tree->Coordinates[dim];
            delta += pow(dx, 2);
        } // for (int dim=0; dim<3; dim++)

        tree->Deltas      = sqrt(delta);
        tree->Coordinates = com;
        tree->Softening   = hmax = &hmax0;
    }

    for(int i = 0; i < 3; i++){
        delete[] quadi[i];
    }
    delete[] quadi;
    delete[] comi;
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

    bool IsParticle = true; // checking whether the current node is a particle
    if ( tree->par == nullptr )   IsParticle = false;

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