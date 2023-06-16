import os
import numpy as np
import matplotlib.pyplot as plt

file_list = [i for i in os.listdir() if "Data_" in i]
Origin      = np.loadtxt('Particle.dat', skiprows=1)
Npar = Origin.shape[0]

if not os.path.exists('./Projection_Par'):
    os.mkdir('./Projection_Par')

for data in file_list:
    image_name = data.split('.')[0]+".png"
    info    = np.loadtxt(data, skiprows=1)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    ax.scatter(info[:, 0], info[:, 1], s=5)
    ax.set_xlim(-1, 2)
    ax.set_ylim(-1, 2)
    plt.savefig('./Projection_Par/'+image_name, dpi=150,bbox_inches='tight')
    plt.close()
