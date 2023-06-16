#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <omp.h>
#include "octree.h"
#include "treewalk.h"
#include "update.h"

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

void WriteAcc(double **g, int Npar, string output)
{
    ofstream out;
    out.open(output);
    out.precision(16);
    out << scientific; // using scientific notation
	for (int i = 0; i < Npar; ++i) {
        if ( i != Npar-1 )
        {
            out << g[i][0] << " ";
            out << g[i][1] << " ";
            out << g[i][2] << endl;
        }
        else
        {
            out << g[i][0] << " ";
            out << g[i][1] << " ";
            out << g[i][2];
        }
        
	}
    out.close();
}

void WriteParticle(Octree *tree, int Npar, double theta, string output)
{
    ofstream out;
    out.open(output);
    out.precision(16);
    out << scientific; // using scientific notation

    out << Npar << " ";
    out << theta << "\n";

	for (int i = 0; i < Npar; ++i) {
        double *pos = tree->partree_arr[i]->par->pos;
        double m    = tree->partree_arr[i]->par->mass;
        double h    = tree->partree_arr[i]->par->softening;

        if ( i != Npar-1 )
        {
            out << pos[0] << " ";
            out << pos[1] << " ";
            out << pos[2] << " ";
            out << m      << " ";
            out << h      << endl;
        }
        else
        {
            out << pos[0] << " ";
            out << pos[1] << " ";
            out << pos[2] << " ";
            out << m      << " ";
            out << h;
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

    double start1 = omp_get_wtime(); 
    Octree* tree = new Octree( Npar, pos, m, h );

    // Declare the array to store the potential
    double* phi = new double[Npar];

    double start2 = omp_get_wtime(); 
    phi = PotentialTarget_tree(Npar, pos, h, tree, G, theta);

    // Declare the array to store the acceleration
    double** g = new double*[Npar];

    double start3 = omp_get_wtime(); 
    g = AccelTarget_tree(Npar, pos, h, tree, G, theta);
    double start4   = omp_get_wtime();


    WritePot(phi, Npar, "./Potential_tree.dat");
    WriteAcc(g, Npar, "./Accel_tree.dat");

    printf("The potential is saved to Potential_tree.dat.\n");
    printf("The acceleration is saved to Accel_tree.dat. \n\n");

    WriteParticle( tree, Npar, theta, "./Timestep00.dat");

    for ( int i = 1; i < 101; i++ ) {
        
        double t = 0.008;
        
        tree_update( tree, t, g );
        for (int j = 0; j < Npar; j++) {
            Particle *cur_par = tree->partree_arr[j]->par;
            for (int k = 0; k < 3; k++)
                pos[j][k] = cur_par->pos[k];
        }
        g = AccelTarget_tree( Npar, pos, h, tree, G, theta );

        if (i % 5 == 0) {
            if (i < 50) {
                WriteParticle( tree, Npar, theta, "./Timestep0" + to_string(i / 5) + ".dat");
                std::cout << " Write  " << i/5 << "th timestep " << t*i << " s file.\n";
            }
            else {
                WriteParticle( tree, Npar, theta, "./Timestep" + to_string(i / 5) + ".dat");
                std::cout << " Write " << i/5 << "th timestep " << t*i << " s file.\n";
            }
        }
    }
    double end   = omp_get_wtime();
    
    printf("\n");
    printf("Wall time for building the tree              = %5.3e s\n", start2 - start1);
    printf("Wall time for calculating the potential      = %5.3e s\n", start3 - start2);
    printf("Wall time for calculating the acceleration   = %5.3e s\n", start4 - start3);
    printf("Wall time for evaluating particles' position = %5.3e s\n", end - start4   );
    printf("\n");
    printf("~ ~ ~ Done ~ ~ ~\n");

    delete tree;
    delete[] m;
    delete[] h;
    delete[] phi;
    for(int i = 0; i < Npar; i++){
        delete[] pos[i];
    }
    delete[] pos;
    for(int i = 0; i < Npar; i++){
        delete[] g[i];
    }
    delete[] g;

    return 0;
}