# Usage: python CreatePar.py Npar
import sys
import numpy as np
from pytreegrav import Potential

N = int(sys.argv[1]) # number of particles
pos = np.random.rand(N,3).astype(np.float64) # positions randomly sampled in the unit cube
m = np.repeat(1./N,N).astype(np.float64) # masses - let the system have unit mass
h = np.repeat(0.01,N).astype(np.float64) # softening radii

f = open("./Particle.dat","w")
f.write("%d\n"%N)

for i in range(N):
    if i==N-1:
        f.write("%.16e %.16e %.16e %.16e %.16e" %(pos[i, 0], pos[i, 1], pos[i, 2], m[i], h[i]))
    else:
        f.write("%.16e %.16e %.16e %.16e %.16e\n" %(pos[i, 0], pos[i, 1], pos[i, 2], m[i], h[i]))

f.close()

# Store reference potential
# Using brute force
p = Potential(pos,m,h,method='bruteforce')
f = open("./Potential_bruteforce_ref.dat", "w")

for i in range(N):
    if i==N-1:
        f.write("%.16e" %(p[i]))
    else:
        f.write("%.16e\n" %(p[i]))

f.close()

# Using tree
p = Potential(pos,m,h,method='tree')
f = open("./Potential_tree_ref.dat", "w")

for i in range(N):
    if i==N-1:
        f.write("%.16e" %(p[i]))
    else:
        f.write("%.16e\n" %(p[i]))

f.close()