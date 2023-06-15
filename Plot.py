import numpy as np
import matplotlib.pyplot as plt

Origin      = np.loadtxt('Particle.dat', skiprows=1)
# Timestep    = np.loadtxt('Timestep.dat', skiprows=1)
Npar = Origin.shape[0]

plt.figure(figsize=(8,6))
plt.scatter(Origin[:, 0], Origin[:, 1])
plt.xlim(-1, 2)
plt.ylim(-1, 2)
plt.savefig("0.png")
# plt.show()

for i in range(1, 10):
    Timestep    = np.loadtxt('Timestep0'+str(i)+'.dat', skiprows=1)
    plt.figure(figsize=(8,6))
    plt.scatter(Timestep[:, 0], Timestep[:, 1])
    plt.xlim(-1, 2)
    plt.ylim(-1, 2)
    plt.savefig(str(i)+".png")
