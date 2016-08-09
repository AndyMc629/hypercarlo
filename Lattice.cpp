#include "Lattice.h"
#include<fstream>
//Lattice constructor
Lattice::Lattice(int a, int b, int c) : Nx(a),Ny(b),Nz(c),lattice(a*b*c) 
{}    
//Lattice destructor
Lattice::~Lattice() {}
//Lattice member func, get volume
int Lattice::Vol(void) { 
	return Nx*Ny*Nz;
}
//Lattice member func, return value at x,y,z coord
Lattice::dipole Lattice::get_dipole(int x,int y,int z) {
	return lattice[x+y*Nx+z*Ny*Nx];
}
//Lattice member func, initialise lattice dep on string
void Lattice::initialise_lattice(std::string s) {
	
	if(s.compare("FERRO")==0) {
	//initialise stuff goes here ....
	std::cout << "FERRO CHOSEN" << std::endl;
	for(int i=0;i<Nx;i++) {
		for(int j=0;j<Ny;j++) {
			for(int k=0;k<Nz;k++) {
				lattice[i+j*Ny+k*Nz*Ny].x = 0;
				lattice[i+j*Ny+k*Nz*Ny].y = 0;
				lattice[i+j*Ny+k*Nz*Ny].z = 1;
				}
			}	
		}
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
//Lattice member func, output lattice to file "string"
void Lattice::output_lattice(std::string datafile) {
std::ofstream output;
output.open(datafile.c_str());
for(int i=0;i<Nx;i++) {
        for(int j=0;j<Ny;j++) {
            for(int k=0;k<Nz;k++) {             
                output << i << "," << j << "," << k << "," << get_dipole(i,j,k).x << ","
                << get_dipole(i,j,k).y << ","
                << get_dipole(i,j,k).z << "\n";
                }
            }
        }
output.close();
}
//Lattice member func, runs simulation on lattice
void Lattice::run() {
/*For equilibriate steps/dipole, choose a site at random and
perform MC step on it*/
/*PSEUDOCODE*/
//for(m=0;m<EquilSteps;m++) {
//(x,y,z)=rand;
//MC_Step(x,y,z);
//}

/*For run steps/dipole, choose a site at random and perform
MC step on it*/
/*PSEUDOCODE*/
//for(n=0;n<RunSteps;m++) {
//(x,y,z)=rand;
//MC_Step(x,y,z);
//}
}
//Lattice member func, performs MC step on x,y,z coordinate dipole
void Lattice::MC_Step(int x, int y, int z) {
/*Recover energy variable, should make that private?
Then propose random new dipole, calc change in energy for that new config
relative to old config, accept or reject new dipole.*/
//Note: see paper http://csml.northwestern.edu/resources/Preprints/mclr.pdf
float E_beforeFlip=site_Hamiltonian(x,y,z);
float theta_new=rand(0,pi),phi_new=rand(0,2*pi);

}
float Lattice::site_Hamiltonian(int x, int y, int z) {
//calculate the energy of a site with coord's x,y,z
	//ISING
	//must calc with three seperate loops so that we don't get
	//next to nearest neighbours from (x+1,y+1,k+1)
	float E;
	float J=1.0;
	for (int i=-1;i<=1;i+=2) {
	E += -J*dot_dipole(get_dipole(x,y,z),get_dipole(x+i,y,z));
	}
	for (int j=-1;j<=1;j+=2) {
	E += -J*dot_dipole(get_dipole(x,y,z),get_dipole(x,y+j,z));		
	}
	for (int k=-1;k<=1;k+=2) {
	E += -J*dot_dipole(get_dipole(x,y,z),get_dipole(x,y,z+k));
	}
return E;
}
float Lattice::dot_dipole(Lattice::dipole p1, Lattice::dipole p2){
float dot=p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
return dot;
}

























