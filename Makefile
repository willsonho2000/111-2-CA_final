all: octree.cpp treewalk.cpp main.cpp
	g++-12 octree.cpp treewalk.cpp update.cpp main.cpp -fopenmp -o main.out
	./main.out Particle.dat

clean:
	rm -f *.out
	rm -f Timestep*.dat
	rm -f *.png
