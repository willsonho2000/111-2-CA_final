import numpy as np
import matplotlib.pyplot as plt

Origin      = np.loadtxt('Particle.dat', skiprows=1)
Npar = Origin.shape[0]

plt.figure(figsize=(8,6))
plt.scatter(Origin[:, 0], Origin[:, 1], s=1)
plt.xlim(-0.5, 1.5)
plt.ylim(-0.5, 1.5)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.text(-0.3, 1.3, "%.3f s" % (0.04 * 0), bbox=props, fontsize=14)

plt.savefig("00.png")
plt.close()

for i in range(1, 21):
    Timestep    = np.loadtxt('Timestep%02d.dat' % (i), skiprows=1)
    plt.figure(figsize=(8,6))
    plt.scatter(Timestep[:, 0], Timestep[:, 1], s=1)
    plt.xlim(-0.5, 1.5)
    plt.ylim(-0.5, 1.5)
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(-0.3, 1.3, "%.3f s" % (0.04 * i), bbox=props, fontsize=14)

    plt.savefig('%02d.png' % (i))
    plt.close()
