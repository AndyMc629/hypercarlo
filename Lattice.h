#ifndef LATTICE_H
#define LATTICE_H
#include<vector>
#include<string>
#include<iostream>

class Lattice {
	public:
		//Initialise a lattice 
		Lattice (int,int,int);
		//Grabs volume of lattice/crystal.
        int Vol(void);
		/*Define dipole struct. NB: next 2 lines
		were in private before but I want user to
		be able to access these, I think.
		i.e means I can change some vectors ad hoc.*/
		struct dipole {float x; float y; float z;};
		std::vector<dipole> lattice;	
		//get the dipole at xyz.
		dipole get_xyz(int,int,int);
		/*initialisation state dependent on string.
		string can have values:
		FERRO = ferro arrangement 1,0,0 or something
		PARA = paraelectric , rand() ...
		PREV = read in previous state as initial state
	
		NB: I'll prob make this private later ...?*/
		void initialise_lattice(std::string);
	
	private:        
		//sizes of crystal
		int Nx,Ny,Nz;
		//void initialise_lattice(std::string);
};

#endif
