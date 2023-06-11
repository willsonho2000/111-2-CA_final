import numpy as np

Pot_tree_ref = np.loadtxt('Potential_tree_ref.dat')
Pot_bru_ref = np.loadtxt('Potential_bruteforce_ref.dat')
Pot_tree = np.loadtxt('Potential_tree.dat')

print("Number of particles = ", Pot_tree.shape[0])
print("\nCalculating error for potential...")
print("Average error between our tree and bruteforce = ", np.mean(np.abs((Pot_bru_ref - Pot_tree)/Pot_bru_ref)))
print("Average error between our tree and pytreegrav = ", np.mean(np.abs((Pot_tree_ref - Pot_tree)/Pot_tree_ref)))
print("Average error between pytreegrav and bruteforce = ", np.mean(np.abs((Pot_tree_ref - Pot_bru_ref)/Pot_bru_ref)))

Accel_tree_ref = np.loadtxt('Accel_tree_ref.dat')
Accel_bru_ref = np.loadtxt('Accel_bruteforce_ref.dat')
Accel_tree = np.loadtxt('Accel_tree.dat')

print("\nCalculating error for acceleration...")
print("Average error between our tree and bruteforce = ", np.mean(np.abs((Accel_bru_ref - Accel_tree)/Accel_bru_ref)))
print("Average error between our tree and pytreegrav = ", np.mean(np.abs((Accel_tree_ref - Accel_tree)/Accel_tree_ref)))
print("Average error between pytreegrav and bruteforce = ", np.mean(np.abs((Accel_tree_ref - Accel_bru_ref)/Accel_bru_ref)))