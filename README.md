# Tree Algorithm For Gravity

## Team Name
Happy Tree Friends

## Collaborators
* Meng-Yuan, Ho (NTU)
* Wei-An, Chen (NTU)

## Presentation
[Unhappy tree](https://docs.google.com/presentation/d/17UjzHAg7b2EwdX5AHwLT_Q2p9SK1uNZrL_Nhc4g7Vlk/edit?usp=sharing)

# Introduction
In this project, we demonstrate a method to calculate the gravitational potential for N-body problem by building up an octree.

# Requirement
* GCC or any C++ compiler
* OpenMP
* Python3 (numpy, pytreegrav)
<p align="right">(<a href="#readme-top">back to top</a>)</p>

# Installation
```
git clone https://github.com/willsonho2000/111-2-CA_final.git
```
<p align="right">(<a href="#readme-top">back to top</a>)</p>

# Implementation

1. First prepare your N particles' information, including position in (x, y, z), mass, and softenting. You can generate a random file using `CreatePar.py`.
```
python3 CreatePar.py Npar theta
```
    where **Npar** is the number of particles and **theta** is the cell-opening criteria. Or, you can manually write a file (_recommand_: Particle.dat) in the following format.
```
Npar theta (first line)
pos_i[x] pos_i[y] pos_i[z] mass_i softening_i
```
    Example:
```
5 6.9999999999999996e-01 
3.2950918101225513e-01 7.6783853139694935e-01 9.7080756233442167e-01 2.0000000000000001e-01 1.0000000000000000e-02
6.4971021323611811e-01 4.9397001999429380e-01 2.8806573895587062e-01 2.0000000000000001e-01 1.0000000000000000e-02
3.2013081665419896e-01 2.3045658135028713e-01 6.8419172781242443e-01 2.0000000000000001e-01 1.0000000000000000e-02
4.1274501434331423e-01 2.4184476920871267e-01 3.0606663461742345e-01 2.0000000000000001e-01 1.0000000000000000e-02
4.9719274486467258e-01 8.8379484672301611e-01 6.3428126058099021e-01 2.0000000000000001e-01 1.0000000000000000e-02
...
```
    _Note: The data should all be saved as **double** type_.

    Potential_tree_ref.dat and Potential_bruteforce_ref.dat will be also generated as reference from [pytreegrav](https://github.com/mikegrudic/pytreegrav).

2. Edit `Makefile` file, where ***mpicxx*** should be replaced by the compiler which the OpenMP is installed in.

3. To compute the potential, you can simply type _make_ to do it.

```
make
```
    or
```
g++-12  octree.cpp treewalk.cpp main.cpp -fopenmp -o main.out
./main.out Particle.dat
```

    Potential_tree.dat will be generated and the potential is stored in the same (particle) order of Particle.dat.

4. Use CalError.py to compute and compare the errors.

```
python3 CalError.py
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

# Performace

<p align="right">(<a href="#readme-top">back to top</a>)</p>
