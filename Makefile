all:
	g++ src/main.cpp src/Lattice.cpp src/mathsFuncs.cpp -std=c++11 -o HyPerCarlo
clean:
	rm HyPerCarlo