#include "Lattice.h"
#include<fstream>
#include<random> //for MT algorithm.
#include "Constants.h" //constants like pi.
#include<cmath>

using namespace Constants; //constants like pi.

std::mt19937 m_rng; //forward declaring m_rng? does this make it work?


//Lattice constructor
Lattice::Lattice(int a, int b, int c, std::mt19937& rng) : Nx(a),Ny(b),Nz(c),m_rng(rng),lattice(a*b*c) 
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
	
	for(int i=0;i<Nx;i++) {
          for(int j=0;j<Ny;j++) {
              for(int k=0;k<Nz;k++) {
                  lattice[i+j*Ny+k*Nz*Ny].x = randomNumber(-1,1);
                  lattice[i+j*Ny+k*Nz*Ny].y = randomNumber(-1,1);
                  lattice[i+j*Ny+k*Nz*Ny].z = randomNumber(-1,1);
                  }
              }
          }

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
void Lattice::Equilibrate(int steps, float T) {
for (int i=0;i<=steps;i++) {
	MC_Step(int(randomNumber(0,Nx)),int(randomNumber(0,Ny)),int(randomNumber(0,Nz)),T); 
}
}
//Lattice member func, runs statistics run on lattice
void Lattice::Run(int steps) {
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
void Lattice::MC_Step(int x, int y, int z, float T) {
/*Recover energy variable, should make that private?
Then propose random new dipole, calc change in energy for that new config
relative to old config, accept or reject new dipole.*/
//Note: see paper http://csml.northwestern.edu/resources/Preprints/mclr.pdf

//Calc current energy and dipole
float E_beforeFlip=site_Hamiltonian(x,y,z);
dipole p_old = lattice[x+y*Nx+z*Nx*Ny];

//generate new trial dipole direction.
float theta = randomNumber(0,pi), phi=randomNumber(0,2*pi); 
dipole p_new;
p_new.x=sin(theta)*cos(phi); //spherical polars.
p_new.y=sin(theta)*sin(phi); //assume r=1 in each case
p_new.z=cos(theta);
//normalise to 1
float norm = sqrt(p_new.x*p_new.x+p_new.y*p_new.y+p_new.z*p_new.z);
p_new.x=p_new.x/norm;p_new.y=p_new.y/norm;p_new.z=p_new.z/norm;
//output for testing
//std::cout << "old (px,py,pz)= " << p_old.x << " " << p_old.y << " " << p_old.z << std::endl;
//std::cout << "new (px,py,pz)= " << p_new.x << " " << p_new.y << " " << p_new.z << std::endl;

//update lattice point as new dipole and calc new energy
lattice[x+y*Nx+z*Nx*Ny]=p_new;
//std::cout << "lattice (px,py,pz)= " << lattice[x+y*Nx+z*Nx*Ny].x << " " << lattice[x+y*Nx+z*Nx*Ny].y << " " << lattice[x+y*Nx+z*Nx*Ny].z << std::endl;
float E_afterFlip=site_Hamiltonian(x,y,z);

//std::cout<< "Ediff = " << E_afterFlip-E_beforeFlip << std::endl;
float dE=E_afterFlip-E_beforeFlip;

//If energy change is positive need Boltzmann factor vs rand number
if (dE>0) {
	float beta=300/(0.025*T);
	float p_boltz=exp( -beta*dE );
	if (p_boltz < randomNumber(0,1)) { //i.e if boltz loses, reject change.
		lattice[x+y*Nx+z*Nx*Ny]=p_old;
	}
}
else if(dE<=0) {
//do nothing, i.e accept the new move.
}
}//end of MC_Step.

float Lattice::site_Hamiltonian(int x, int y, int z) {
//calculate the energy of a site with coord's x,y,z
	//ISING
	//must calc with three seperate loops so that we don't get
	//next to nearest neighbours from (x+1,y+1,k+1)
	float E;
	float J=0.025;
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
//create random number between min and max using mersenneTwister.
float Lattice::randomNumber(float min, float max){
        std::uniform_real_distribution<float> uniformDistribution(min, max);
        return uniformDistribution(m_rng);
    }





















