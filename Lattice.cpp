#include "Lattice.h"

//Lattice constructor
Lattice::Lattice(int a, int b, int c) : Nx(a),Ny(b),Nz(c),lattice(a*b*c) 
{}    
//Lattice member func, get volume
int Lattice::Vol(void) { 
	return Nx*Ny*Nz;
}
//Lattice member func, return value at x,y,z coord
Lattice::dipole Lattice::get_xyz(int x,int y,int z) {
	return lattice[x+y*Nx+z*Ny*Nx];
}
//Lattice member func, initialise lattice dep on string
void Lattice::initialise_lattice(std::string s) {
	
	if(s.compare("FERRO")==0) {
	//initialise stuff goes here ....
	std::cout << "FERRO CHOSEN" << std:endl;
	}
	else if(s.compare("PARA")==0) {
	//initialise stuff goes here ....
	std::cout << "PARA CHOSEN" << std::endl;	
	}
	else if(s.compare("PREV")==0) {//previous output
	//initialise stuff goes here ....
	std::cout << "PREV CHOSEN" << std::endl;
	}
}

