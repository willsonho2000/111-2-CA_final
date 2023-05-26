all: octree.cpp treewalk.cpp main.cpp
	g++-12 octree.cpp treewalk.cpp main.cpp -o main.out
	./main.out Particle.dat

clean:
	rm -f *.out
