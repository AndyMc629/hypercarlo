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
int Nx=5,Ny=5,Nz=5;
std::mt19937 rng{std::chrono::high_resolution_clock::now().time_since_epoch().count()};
Lattice lattice=Lattice(Nx,Ny,Nz,rng);
//check it's volume
int vol=lattice.Vol();
/*this initialisation step will eventually happen internally maybe? 
Maybe not I quite like it this way ....*/
lattice.initialise_lattice("FERRO");
lattice.output_lattice("InitialState.dat");

float answer = lattice.dot_dipole(lattice.get_dipole(1,1,1),lattice.get_dipole(1,2,1));
std::cout << answer << std::endl;

for(int i=0;i<10;i++){

lattice.MC_Step(1,1,1);

}

/*for(int i=0;i<Nx;i++) {
        for(int j=0;j<Ny;j++) {
            for(int k=0;k<Nz;k++) {				
                std::cout << lattice.get_xyz(i,j,k).x << " "
               	<< lattice.get_xyz(i,j,k).y << " "
				<< lattice.get_xyz(i,j,k).z << std::endl;
                //std::cout << lattice[i+j*Ny+k*Nz*Ny].z << " ";
                }
            }
        }
*/    
return 0;
}
