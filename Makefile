all:
	g++ src/main.cpp src/Lattice.cpp src/mathsFuncs.cpp src/inputParser.cpp -std=c++11 -o HyPerCarlo
clean:
	rm HyPerCarlo *.dat
