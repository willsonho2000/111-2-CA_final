#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <algorithm>
using namespace std;

const int N = 5; // number of particles

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

double PotentialWalk_quad(double* pos, tree tree, double softening, double theta)
{
    int no = tree.NumParticles;
    double phi = 0, dx[3] = {};
    double quad[3][3], r5inv;

    while ( no > -1 )
    {
        double r=0;
        for (int k=0; k<3; k++)
        {
            dx[k] = tree.Coordinates[no][k] - pos[k];
            r+=dx[k]*dx[k];
        }
        r = sqrt(r);
        double h = fmax(tree.Softenings[no], softening);
        
        if ( no < tree.NumParticles )
        {
            if ( r > 0 )
            {
                if ( r < h)
                {
                    phi += tree.Masses[no] * PotentialKernel(r,h); 
                }
                else
                {
                    phi -= tree.Masses[no] / r;
                } // if ( r < h)
            } // if ( r > 0 )
            no = tree.NextBranch[no];
        } // if ( no < tree.NumParticles )
        else if ( r > fmax(tree.Sizes[no]/theta + tree.Deltas[no], h+tree.Sizes[no]*0.6+tree.Deltas[no]) )
        {
            phi -= tree.Masses[no]/r;
            copy(&tree.Quadrupoles[no][0][0], &tree.Quadrupoles[no][0][0] + 3 * 3, &quad[0][0]);
            r5inv = 1 / pow(r, 5);

            for (int k=0; k<3; k++)
            for (int l=0; l<3; l++)
            phi -= 0.5 * dx[k] * quad[k][l] * dx[l] * r5inv;

            no = tree.NextBranch[no];
        } // else if ( r > fmax(tree.Sizes[no]/theta + tree.Deltas[no], h+tree.Sizes[no]*0.6+tree.Deltas[no]) )
        else
        {
            no = tree.FirstSubnode[no];
        } // else
    } // while ( no > -1 )

    return phi;
}

double* PotentialTarget_tree(double** pos_target, double* softening_target, tree tree, int G, double theta)
{
    double result[N];

    for (int i=0; i<N; i++)
    result[i] = G*PotentialWalk_quad(pos_target[i], tree, softening_target[i], theta);

    return result;
}


int main(){
    int *arr = fun();
    cout << arr[0] << endl;


    return 0;
}