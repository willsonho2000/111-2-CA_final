all: octree.cpp treewalk.cpp main.cpp
	mpicxx  octree.cpp treewalk.cpp main.cpp -o main.out
	./main.out Particle.dat

clean:
	rm -f *.out
