# Tree Algorithm For Gravity

## Team Name
Happy Tree Friends

## Collaborators
* Meng-Yuan Ho (NTU)
* Wei-An Chen (NTU)

## Presentation
[Happy tree](https://docs.google.com/presentation/d/17UjzHAg7b2EwdX5AHwLT_Q2p9SK1uNZrL_Nhc4g7Vlk/edit?usp=sharing)

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
where **Npar** is the number of particles and **theta** is the cell-opening criteria. Or, you can manually write a file (_recommend_: `Particle.dat`) in the following format.
```
Npar theta (first line)
pos_i[x] pos_i[y] pos_i[z] vel[x] vel[y] vel[z] mass_i softening_i
```
Example:
```
5 6.9999999999999996e-01 
8.3418945034164693e-01 3.4578865673568515e-01 5.3318914437914633e-01 -1.4892456335172388e-01 1.8642446499586351e-01 1.7642995050771193e-01 1.0000000000000000e-02 1.0000000000000000e-02
2.5739577825541171e-01 3.8865321100378047e-01 9.6849050947488047e-01 -8.9216643422320718e-03 1.5872282261062354e-01 -3.5179269862034235e-02 1.0000000000000000e-02 1.0000000000000000e-02
8.5648279277900918e-01 4.8043861303190716e-01 6.1264434494372688e-01 2.1333913747113697e-01 -1.6720246137579009e-01 -6.8057000867475403e-02 1.0000000000000000e-02 1.0000000000000000e-02
5.5746719784732779e-02 9.7598227534745652e-01 2.5970268980899947e-01 1.5212239525997318e-01 -1.0677190282585058e-01 1.5485664340995503e-01 1.0000000000000000e-02 1.0000000000000000e-02
5.4413433587357307e-01 1.5550535592805281e-01 2.6869243960931455e-01 -2.4509868726131129e-01 -1.7140646768390189e-01 -1.7696595983169527e-01 1.0000000000000000e-02 1.0000000000000000e-02
...
```
_Note: The data should all be saved as **double** type_.

`Potential_tree_ref.dat`, `Potential_bruteforce_ref.dat`, `Accel_bruteforce_ref.dat` and `Accel_tree_ref.dat` will be also generated as reference from [pytreegrav](https://github.com/mikegrudic/pytreegrav).

2. Edit `Makefile`, where ***mpicxx*** should be replaced by the compiler which the OpenMP is installed in.

3. To compute the potential, you can simply type _make_ to do it.

```
make
```
or
```
g++-12  octree.cpp treewalk.cpp main.cpp -fopenmp -o main.out
./main.out Particle.dat
```

`Potential_tree.dat` and `Accel_tree.dat` will be generated and the order is stored in the same (particle) as `Particle.dat`.

4. Use `CalError.py` to compute and compare the errors.

```
python3 CalError.py
```

5. Use `Plot.py` to plot projection and `movie.sh` to make movie.

```
python3 Plot.py
sh movie.sh
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

