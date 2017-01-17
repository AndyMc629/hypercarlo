/* Source code for model class, will contain the Hamiltonians for different 
 * models etc.
 */
#include "Model.h"

//Model constructor
Model::Model(std::string chosenModel) : modelName(chosenModel) {
    std::cout << "Selected model: " << chosenModel << " model\n";
}
//Model destructor
Model::~Model() {}

double Model::Ising_Site_Hamiltonian(int x, int y) {
    //ISING
    //must calc with three seperate loops so that we don't get
    //next to nearest neighbours from (x+1,y+1,k+1)
    double H=0.0;
    for (int i=-1;i<=1;i+=2) {
    H += -J*Lattice::dot_dipole(get_dipole(x,y),get_dipole(x+i,y));
    }
    for (int j=-1;j<=1;j+=2) {
    H += -J*Lattice::dot_dipole(get_dipole(x,y),get_dipole(x,y+j));		
    }
    return H;    
}