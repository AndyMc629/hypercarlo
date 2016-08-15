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
//Note that for PBC have to find mod(Ni-1) for each
//coordinate (minus one
//because c++ arrays start at zero).
Lattice::dipole Lattice::get_dipole(int x,int y,int z) {
	return lattice[(x%(Nx-1))+(y%(Ny-1))*Nx+(z%(Nz-1))*Ny*Nx];
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
	
    else if(s.compare("PARA_ISING")==0) {
	//initialise stuff goes here ....
	std::cout << "PARA ISING CHOSEN" << std::endl;	
	float random=randomNumber(0,1);
	for(int i=0;i<Nx;i++) {
          for(int j=0;j<Ny;j++) {
              for(int k=0;k<Nz;k++) {
                  lattice[i+j*Ny+k*Nz*Ny].x = 0; 
                  lattice[i+j*Ny+k*Nz*Ny].y = 0; 
                  random=randomNumber(0,1);
				if(random>=0.5){
					lattice[i+j*Ny+k*Nz*Ny].z = 1;
				}
                  else{
                  lattice[i+j*Ny+k*Nz*Ny].z = -1;
                  }    
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
                output << i << " " << j << " " << k << " " << get_dipole(i,j,k).x << " "
                << get_dipole(i,j,k).y << " "
                << get_dipole(i,j,k).z << "\n";
                }
            }
        }
output.close();
}
void Lattice::Equilibrate(int stepsPerSite, float T) {
//Open data file to store equilibration stats.
std::ofstream output;
output.open("EquilibrationStats.dat");
output << "#Equilibrated stats.\n"
		<< "#kT/J nstep/site E_tot/site P_tot/site\n";   
for (int i=0;i<=(stepsPerSite*Vol());i++) {	
	MC_Step(int(randomNumber(0,Nx)),int(randomNumber(0,Ny)),int(randomNumber(0,Nz)),T); 
	if(((10*i)%(stepsPerSite*Vol()))==0){ //every 10% output something
		output << (0.025*T)/(J*(float)300) << " " <<  (float)i/(Vol()) << " " 
		<< total_Energy()/(Vol()) << " " << total_Polarisation()/(Vol()) << "\n"; 
		std::cout << ((float)i/stepsPerSite) << " steps/site: E/vol = " 
		<< total_Energy()/Vol() << ", P/vol = " << total_Polarisation()/Vol() << std::endl;
	}	
}
output.close();
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
dipole p_new;
p_new.x=0; //Ising so only z component.
p_new.y=0; //Ising so only z component.
p_new.z= -p_old.z; //Ising so can only flip.
//normalise to 1 - no need in Ising.
//float norm = sqrt(p_new.x*p_new.x+p_new.y*p_new.y+p_new.z*p_new.z);
//p_new.x=p_new.x/norm;p_new.y=p_new.y/norm;p_new.z=p_new.z/norm;

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
	float E=0.0;
	//float J=0.025;
	for (int i=-1;i<=1;i+=2) {
	E += -J*dot_dipole(get_dipole(x,y,z),get_dipole(x+i,y,z));
	}
	for (int j=-1;j<=1;j+=2) {
	E += -J*dot_dipole(get_dipole(x,y,z),get_dipole(x,y+j,z));		
	}
	for (int k=-1;k<=1;k+=2) {
	E += -J*dot_dipole(get_dipole(x,y,z),get_dipole(x,y,z+k));
	}
if(E>10e2) {
std::cout << "\n \n \n ********************************\n"
			<<" WARNING! High energy on a site, most likely an error. More info below:" 
			<<" \n In site H func: \n"; 
for (int i=-1;i<=1;i+=2) {
     E += -J*dot_dipole(get_dipole(x,y,z),get_dipole(x+i,y,z));
     std::cout << "dot1= "<< -J*dot_dipole(get_dipole(x,y,z),get_dipole(x+i,y,z)) << "\n";
		}
     for (int j=-1;j<=1;j+=2) {
     E += -J*dot_dipole(get_dipole(x,y,z),get_dipole(x,y+j,z));
	std::cout << "dot2 = "<< -J*dot_dipole(get_dipole(x,y,z),get_dipole(x,y+j,z)) << "\n";
     }
     for (int k=-1;k<=1;k+=2) {
     E += -J*dot_dipole(get_dipole(x,y,z),get_dipole(x,y,z+k));
	std::cout << "dot3 = "<< -J*dot_dipole(get_dipole(x,y,z),get_dipole(x,y,z+k)) << "\n";
     }

}
return E;
}
float Lattice::total_Energy(){
float E=0.0;
float E_site=0.0;
 for(int i=0;i<Nx;i++) {
          for(int j=0;j<Ny;j++) {
              for(int k=0;k<Nz;k++) {
                 	E_site=site_Hamiltonian(i,j,k); 
					E += E_site;
				if(E_site>10e2) {
				std::cout << "\n \n \n ********************************\n" 
		             <<" WARNING! High energy on a site, most likely an error. Lattice config will be outputted\n"
						<<"More info below:\n";
					std::cout << "E>100 = " << site_Hamiltonian(i,j,k) <<" at site ("<< i << ", "<< j << ", "<< k <<")" << std::endl;
					output_lattice("ProblemLattice.dat");
					}
				}
              }
          }      
return E;
}
//Calc. total polarisation of lattice.
//Choice of words not sacrosanct, this also
//refers to the magnetisation for spin simulations.
float Lattice::total_Polarisation(){
//vector to store net pol.
dipole P_net;P_net.x=0.0;P_net.y=0.0;P_net.z=0.0;
	for(int i=0;i<Nx;i++) {
    	      for(int j=0;j<Ny;j++) {
        	      for(int k=0;k<Nz;k++) {
					P_net.x += lattice[i+j*Nx+k*Nx*Ny].x;
					P_net.y += lattice[i+j*Nx+k*Nx*Ny].y;
					P_net.z += lattice[i+j*Nx+k*Nx*Ny].z;	
				}
			}
		}
float P_net_mag=sqrt(dot_dipole(P_net,P_net));
return P_net_mag;
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





















