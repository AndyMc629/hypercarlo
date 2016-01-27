/* Monte Carlo lattice model simulator, this will be generic and 3D. All that will need changed/updated in future will be Hamiltonians!

Author: Andrew P. McMahon, Imperial College London.
Began: 25/1/2016
*/

#include<cstdlib>
#include<iostream>
#include "Lattice.h"
#include<vector>

using namespace std;

int main() {
//Set size of lattice
Lattice lattice=Lattice(10,10,10);
//check it's volume
int vol=lattice.Vol();
/*this initialisation step will eventually happen internally maybe? 
Maybe not I quite like it this way ....*/
lattice.initialise_lattice("PARA");

/*Tests for accessinv volume and lattice data*/
//cout << "Volume of lattice: " << vol << endl;
//Lattice::dipole p1 = lattice.get_xyz(1,2,3);
//cout << "L[1][2][3]: " << p1.x << endl;


return 0;
}
