#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <algorithm>

void ComputeMoments(double* mass, double com[3], double* hmax, Octree* tree);

double PotentialKernel(double r, double h);

double PotentialWalk_quad(double* pos, Octree* tree, double theta, double softening);
double* PotentialTarget_tree(int Npar, double** pos_target, double* softening_target, Octree* tree, int G, double theta);
