all: octree.cpp test.cpp
	g++ octree.cpp test.cpp -o test.out
	./test.out

clean:
	rm -f *.out
