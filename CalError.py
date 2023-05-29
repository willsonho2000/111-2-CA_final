import numpy as np

tree_ref = np.loadtxt('Potential_tree_ref.dat')
tree_bru = np.loadtxt('Potential_bruteforce_ref.dat')
tree = np.loadtxt('Potential_tree.dat')

print("Number of particles = ", tree_ref.shape[0])
print("Average error between our tree and bruteforce = ", np.mean(np.abs((tree_bru - tree)/tree_bru)))
print("Average error between our tree and pytreegrav = ", np.mean(np.abs((tree_ref - tree)/tree_ref)))
print("Average error between pytreegrav and bruteforce = ", np.mean(np.abs((tree_ref - tree_bru)/tree_bru)))