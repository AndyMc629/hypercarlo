/* Monte Carlo lattice model simulator, this will be generic and 3D. All that will need changed/updated in future will be Hamiltonians!

Author: Andrew P. McMahon, Imperial College London.
Began: 25/1/2016
*/

#include<cstdlib>
#include<iostream>
#include "Lattice.h"
#include<vector>

//using namespace std;

int main() {
//Set size of lattice
int Nx=5,Ny=5,Nz=5;
Lattice lattice=Lattice(Nx,Ny,Nz);
//check it's volume
int vol=lattice.Vol();
/*this initialisation step will eventually happen internally maybe? 
Maybe not I quite like it this way ....*/
lattice.initialise_lattice("FERRO");
lattice.output_lattice("InitialState.dat");

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
