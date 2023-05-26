# Usage: python CreatePar.py Npar
import sys
import numpy as np
from pytreegrav import Potential

N = int(sys.argv[1]) # number of particles
pos = np.random.rand(N,3) # positions randomly sampled in the unit cube
m = np.repeat(1./N,N) # masses - let the system have unit mass
h = np.repeat(0.01,N) # softening radii

f = open("./Particle.dat","w")

for i in range(N):
    if i==N-1:
        f.write("%13.7e %13.7e %13.7e %13.7e %13.7e" %(pos[i, 0], pos[i, 1], pos[i, 2], m[i], h[i]))
    else:
        f.write("%13.7e %13.7e %13.7e %13.7e %13.7e\n" %(pos[i, 0], pos[i, 1], pos[i, 2], m[i], h[i]))

f.close()

# Store reference potential
p = Potential(pos,m,h)
f = open("./Ref.dat", "w")

for i in range(N):
    if i==N-1:
        f.write("%13.7e" %(p[i]))
    else:
        f.write("%13.7e\n" %(p[i]))

f.close()