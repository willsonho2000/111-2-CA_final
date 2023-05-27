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
* Python3
<p align="right">(<a href="#readme-top">back to top</a>)</p>

# Installation
```
git clone https://github.com/willsonho2000/111-2-CA_final.git
```
<p align="right">(<a href="#readme-top">back to top</a>)</p>

# Implementation

1. First prepare your N particles' information, including position in (x, y, z), mass, and softenting. You can generate a random file using `CreatePar.py`.
```
python3 CreatePar.py Npar
```
where **Npar** is the number of particles. Or, you can manually write a file (_recommand_: Particle.dat) in the following format.
```
pos_i[x] pos_i[y] pos_i[z] mass_i softening_i
```
_Note: The data should all be saved as **double** type_.

2. Edit `Makefile` file, where ***mpicxx*** should be replaced by the compiler which the OpenMP is installed in.

3. Edit `main.cpp`. Change the Npar to the **exact** number of particle for your problem.

```
const int Npar = 10;    // the number of particles
```

4. To compute the potential, you can simply type _make_ to do it.

```
make
```
or
```
./main.out Particle.dat
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>
