#ifndef LATTICE_H
#define LATTICE_H
#include<vector>
#include<string>
#include<iostream>
#include<random>

class Lattice {
	public:
		//Initialise a lattice/lattice constructor 
		Lattice (int,int,int);
		//Lattice destructor
		~Lattice();
		//Grabs volume of lattice/crystal.
	        int vol(void);
		/*Define dipole struct. NB: next 2 lines
		were in private before but I want user to
		be able to access these, I think.
		i.e means I can change some vectors ad hoc.*/
		struct dipole {float x; float y; float z;};
		//the lattice object is a vector of dipole structs.
		std::vector<dipole> lattice;	
		//get the dipole at xyz.
		dipole get_dipole(int,int,int);
		/*initialisation state dependent on string.
		string can have values:
		FERRO = ferro arrangement 1,0,0 or something
		PARA = paraelectric , rand() ...
		PREV = read in previous state as initial state
		NB: I'll prob make this private later ...?*/
		void initialise_lattice(std::string);
		/*Output to file supplied by string*/
		void output_lattice(std::string);
		/*Equilibrate lattice for int steps/dipole*/
		void equilibrate(int);
		/*Run sim gather statistics for int steps/dipole*/
		void run(int);
	private:        
		//sizes of crystal
		int Nx,Ny,Nz;
		//void initialise_lattice(std::string);
		void mc_step(int,int,int);
		//init a random number for use later ...
		std::mt19937 rng(std::random_device{}());	
};

#endif
