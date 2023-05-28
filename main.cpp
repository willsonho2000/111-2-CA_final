#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include "octree.h"
#include "treewalk.h"

using namespace std;

void ReadPar(int Npar, double** pos, double* m, double* h, string input)
{
    ifstream in;
    in.open(input);
    if(in.fail()){ 
        cout << "input file opening failed";
        exit(1);
    }

    int dum;
    double dumf; // skip Npar and thetae
    in >> dum;
    in >> dumf;

    for (int i=0; i<Npar; i++)
    {
        in >> pos[i][0];
        in >> pos[i][1];
        in >> pos[i][2];
        in >> m[i];
        in >> h[i];
    }

    in.close();
    return;
}

void WritePot(double *phi, int Npar, string output)
{
    ofstream out;
    out.open(output);
    out.precision(16);
    out << scientific; // using scientific notation
	for (int i = 0; i < Npar; ++i) {
        if ( i != Npar-1 )
        {
            out << phi[i] << endl;
        }
        else
        {
            out << phi[i];
        }
        
	}
    out.close();
}

int main( int argc, char* argv[] ) {
    // Use ./main.out ./Particle.dat
    // Basic settings
    string input = argv[1];
    const int    G     =   1;
    const int    NThread = 4;
    omp_set_num_threads( NThread );

    // Read the number of particle
    int Npar;
    double theta;

    ifstream in;
    in.open(input);
    in >> Npar;
    in >> theta;
    in.close();

    printf("Number of particles (Npar)    = %d\n", Npar );
    printf("Cell-opening criteria (theta) = %f\n", theta);
    printf("\n");

    // Declare particle's properties
    double** pos = new double*[Npar];
    for (int i=0; i<Npar; i++) pos[i] = new double[3];

    double* m = new double[Npar];
    double* h = new double[Npar];

    // Store particles' properties
    ReadPar(Npar, pos, m, h, input);

    Octree* tree = new Octree( Npar, pos, m, h );

    // Declare the array to store the potential
    double* phi = new double[Npar];

    double start = omp_get_wtime();
    phi = PotentialTarget_tree(Npar, pos, h, tree, G, theta);
    double end   = omp_get_wtime();

    printf("Wall time = %5.3e s\n", end - start);
    printf("\n");

    WritePot(phi, Npar, "./Potential_tree.dat");

    printf("The potential is saved to Potential_tree.dat.\n");
    printf("~ ~ ~ Done ~ ~ ~\n");

    delete tree;
    delete[] m;
    delete[] h;
    delete[] phi;
    for(int i = 0; i < Npar; i++){
        delete[] pos[i];
    }
    delete[] pos;

    return 0;
}