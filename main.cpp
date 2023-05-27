#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include "octree.h"
#include "treewalk.h"

using namespace std;

void ReadPar(int Npar,double** pos, double* m, double* h, string input)
{
    ifstream in;
    in.open(input);
    if(in.fail()){ 
        cout << "input file opening failed";
        exit(1);
    }

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

int main( int argc, char* argv[] ) {
    // Use ./main.out ./Particle.dat
    // Basic settings
    string input = argv[1];
    const int Npar = 10;

    // Declare particle's properties
    double** pos = new double*[Npar];
    for (int i=0; i<Npar; i++) pos[i] = new double[3];

    double* m = new double[Npar];
    double* h = new double[Npar];

    // Store particles' properties
    ReadPar(Npar,pos, m, h, input);

    double total_mass = 0.;
    for ( int i = 0; i < Npar; i++ ) total_mass += m[i];

    Octree* tree = new Octree( Npar, pos, m, h );
    cout << "total masses: " << tree->Masses << " , compared to the initial condition: " << total_mass << "\n";

    return 0;
}