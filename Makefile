all: octree.cpp treewalk.cpp main.cpp
	mpicxx octree.cpp treewalk.cpp update.cpp main.cpp -fopenmp -o main.out
	./main.out Particle.dat

clean:
	rm -f *.out
