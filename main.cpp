/* Monte Carlo lattice model simulator, this will be generic and 3D. All that will need changed/updated in future will be Hamiltonians!

Author: Andrew P. McMahon, Imperial College London.
Began: 25/1/2016

I imagine the workflow of the program:

Initialise a lattice 
|--> Choose model Hamiltonian? Should I create a "Model" class containing the Hamiltonian and Monte Carlo steps 
     or do I just make this a method in the lattice class?
V
Run equilibration on lattice --> 
|
V
Run simulation on lattice, gather statistics --> Calculate averages --> Output averages.


*/

#include<cstdlib>
#include<iostream>
#include "Lattice.h"
#include<vector>
#include<random>
#include<chrono>
//using namespace std;

int main() {
//Set size of lattice
int Nx=10,Ny=10,Nz=10;
//Seed a Mersenne Twister random number generator.
std::mt19937 rng{std::chrono::high_resolution_clock::now().time_since_epoch().count()};
//Initialise a lattice object.
Lattice lattice=Lattice(Nx,Ny,Nz,rng);


//Initialise the lattice to FERRO or PARA
lattice.initialise_lattice("PARA_ISING");
//output initial lattice.
lattice.output_lattice("InitialState.dat");

//equilibrate lattice at temp T.
float T=1200; //K
int equilStepsPerSite=1000;
lattice.Equilibrate(equilStepsPerSite,T);

//output equilibrated lattice.
//TEST
std::string equilFile="Lattice_equil_" + std::to_string(equilStepsPerSite)+"_stepsPerSite.dat";
lattice.output_lattice(equilFile);

return 0;
}
