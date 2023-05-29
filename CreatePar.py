# Usage: python CreatePar.py Npar theta
import sys
import numpy as np
from pytreegrav import Potential

def NoNpar(arg):
    if len(arg) <= 2:
        raise IndexError("You need to specify the values of Npar and theta")

def NotCorrectVal(arg):
    if int(arg[1]) <= 1:
        raise ValueError("Please enter a correct value of Npar (Npar >= 2)")
    elif float(arg[2]) < 0:
        raise ValueError("Please enter a correct value of theta (theta > 0)")

try:
    NoNpar(sys.argv)
except IndexError as error:
    print(error)
    sys.exit()

try:
    NotCorrectVal(sys.argv)
except ValueError as error:
    print(error)
    sys.exit()

N = int(sys.argv[1]) # number of particles
theta = float(sys.argv[2]) # number of particles

print("Npar = ", N)
print("theta = ", theta)

pos = np.random.rand(N,3).astype(np.float64) # positions randomly sampled in the unit cube
m = np.repeat(1./N,N).astype(np.float64) # masses - let the system have unit mass
h = np.repeat(0.01,N).astype(np.float64) # softening radii

print(" Generating Particle.dat... ")
f = open("./Particle.dat","w")
f.write("%d %.16e \n"%(N, theta))

for i in range(N):
    if i==N-1:
        f.write("%.16e %.16e %.16e %.16e %.16e" %(pos[i, 0], pos[i, 1], pos[i, 2], m[i], h[i]))
    else:
        f.write("%.16e %.16e %.16e %.16e %.16e\n" %(pos[i, 0], pos[i, 1], pos[i, 2], m[i], h[i]))

f.close()
print("Particle.dat is saved.\n")

# Store reference potential
# Using brute force
print("Calculating the potential with brute force...")
p = Potential(pos,m,h,method='bruteforce')
f = open("./Potential_bruteforce_ref.dat", "w")

for i in range(N):
    if i==N-1:
        f.write("%.16e" %(p[i]))
    else:
        f.write("%.16e\n" %(p[i]))

f.close()
print("Potential_bruteforce_ref.dat is saved.\n")

# Using tree
print("Calculating the potential with pytreegrav.")
p = Potential(pos,m,h,method='tree')
f = open("./Potential_tree_ref.dat", "w")

for i in range(N):
    if i==N-1:
        f.write("%.16e" %(p[i]))
    else:
        f.write("%.16e\n" %(p[i]))

f.close()
print("Potential_tree_ref.dat is saved.\n")