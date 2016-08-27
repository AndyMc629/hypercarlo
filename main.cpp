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
#include<fstream>
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
float temp; //K
int equilStepsPerSite=10000;
int ensemble_size=1000;


std::ofstream mainOutput;
mainOutput.open("Output.dat");

mainOutput << "#T(K) E_av Esqrd_av P_av Psqrd_av Cv\n"; 

for (int T=900;T>=50;T=T-50) {
temp=(float)T;
lattice.Equilibrate(equilStepsPerSite,temp);
//output equilibrated lattice.
std::string equilFile="Lattice_equil_" + std::to_string(equilStepsPerSite)+"_stepsPerSite_"+std::to_string((int)T)+"K.dat";
lattice.output_lattice(equilFile);
//calc estimators
lattice.Run(ensemble_size,temp);

std::cout << "T="<<T<<"K:\n"
<< "E_av="<<lattice.E_av << "\n"
<< "Esqrd_av="<<lattice.Esqrd_av << "\n"
<< "Cv="<<lattice.Cv << "\n"
<< "P_av="<<lattice.P_av << "\n"
<< "Psqrd_av="<<lattice.Psqrd_av<<"\n\n";

mainOutput << T << " " << lattice.E_av << " " << lattice.Esqrd_av << " "
<<lattice.P_av<< " " << lattice.Psqrd_av << " " << lattice.Cv << "\n";

}
mainOutput.close();
return 0;
}
