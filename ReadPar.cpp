#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
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
}

int main(int argc, char* argv[]){
    // Use ./ReadPar ./Particle.dat
    string input = argv[1];
    const int Npar = 10;

    double** pos = new double*[Npar];
    for (int i=0; i<Npar; i++) pos[i] = new double[3];

    double* m = new double[Npar];
    double* h = new double[Npar];

    ReadPar(Npar,pos, m, h, input);
    for (int i=0; i<Npar; i++)
    {
        printf("%13.7e %13.7e %13.7e %13.7e %13.7e\n", pos[i][0], pos[i][1], pos[i][2], m[i], h[i]);
    }

    delete[] m;
    delete[] h;
    for(int i = 0; i < Npar; i++){
        delete[] pos[i];
    }
    delete[] pos;

    return 0;
}