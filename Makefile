all: octree.cpp treewalk.cpp test.cpp
	g++-12 octree.cpp treewalk.cpp test.cpp -o test.out
	./test.out

clean:
	rm -f *.out
