/* Source code for the 'Lattice' class, which contains main bulk of run 
 * functions.
 */

#include "Lattice.h"
#include<fstream>
#include<random> //for MT algorithm.
#include "Constants.h" //constants like pi,kB etc.
#include<cmath>
#include "mathsFuncs.h" //for correct modulo function etc...
using namespace Constants; //constants like pi.

//Seed a Mersenne Twister random number generator with current time.
//std::mt19937 m_rng; //forward declaring m_rng? does this make it work?

//Lattice constructor
//Lattice::Lattice(int a, int b, std::mt19937& rng) : Nx(a),Ny(b),m_rng(rng),lattice(a*b) 
Lattice::Lattice(int a, int b, std::string modelName) : Nx(a),Ny(b),model(modelName),lattice(a*b)
{ 
    std::cout << "A lattice object of size "<<Nx<<" x "<<Ny
        <<" has been initialised\n";
}    
//Lattice default constructor
Lattice::Lattice() {}
//Lattice destructor
Lattice::~Lattice() {}
//Lattice member func, get volume
int Lattice::Vol(void) { 
	return Nx*Ny;
}
//Lattice member func, return value at x,y,z coord
//Note that for PBC have to find mod(Nx) for each
//coordinate. Have double checked that it is Nx and not Nx-1.
Lattice::dipole Lattice::get_dipole(int x,int y) {
        return lattice[mod(x,Nx)+mod(y,Ny)*Nx]; //periodic boundary conditions
}
//Lattice member func, initialise lattice dep on string
void Lattice::initialise_lattice(std::string s) {
	
	if(s.compare("FERRO")==0) {
	//initialise stuff goes here ....
	std::cout << "FERRO CHOSEN" << std::endl;
	for(int i=0;i<Nx;i++) {
		for(int j=0;j<Ny;j++) {
				lattice[i+j*Nx].x = 0;
				lattice[i+j*Nx].y = 0;
				lattice[i+j*Nx].z = 1;
			}	
		}
	}
	else if(s.compare("PARA")==0) {
	//initialise stuff goes here ....
	std::cout << "PARA CHOSEN" << std::endl;	
	
	for(int i=0;i<Nx;i++) {
          for(int j=0;j<Ny;j++) {
                  lattice[i+j*Nx].x = randomNumber(-1,1);
                  lattice[i+j*Nx].y = randomNumber(-1,1);
                  lattice[i+j*Nx].z = randomNumber(-1,1);
              }
          }
	}
	
    else if(s.compare("PARA_ISING")==0) {
	//initialise stuff goes here ....
	std::cout << "PARA ISING CHOSEN" << std::endl;	
	double random=randomNumber(0,1);
	for(int i=0;i<Nx;i++) {
          for(int j=0;j<Ny;j++) {
                  lattice[i+j*Nx].x = 0; 
                  lattice[i+j*Nx].y = 0; 
                  random=randomNumber(0,1);
				if(random>=0.5){
					lattice[i+j*Nx].z = 1;
				}
                  else{
                  lattice[i+j*Nx].z = -1;
                	}
                  }
              }
          }	

    else if(s.compare("PREV")==0) {//previous output
	//read in previous output file stuff goes here ....
	std::cout << "PREV CHOSEN" << std::endl;
	}
}
//Lattice member func, output lattice to file "string"
void Lattice::output_lattice(std::string datafile) {
std::ofstream output;
output.open(datafile.c_str());
for(int i=0;i<Nx;i++) {
        for(int j=0;j<Ny;j++) {
                output << i << " " << j << " " << get_dipole(i,j).x << " "
                << get_dipole(i,j).y << " "
                << get_dipole(i,j).z << "\n";
                }
            }
output.close();
}

//Lattice member func, equilibrate for stepsPerSite MCSweeps.
void Lattice::Equilibrate(int stepsPerSite, double T) {
E_total=total_Energy();
P_total=total_Polarisation();
std::cout << "Initial lattice E = " << E_total
        << "\n Initial lattice P = " << P_total << std::endl;
//Open data file to store equilibration stats.
std::ofstream output;
output.open("EquilibrationStats_"+std::to_string(T)+"K.dat");
output << "#Equilibrated stats.\n"
		<< "#kT/J MCS E_tot/site E_tot(global)/site P_tot/site\n";   

for (int i=0;i<=(stepsPerSite*Vol());i++) {	
        MC_Step(int(randomNumber(0,Nx)),int(randomNumber(0,Ny)),T); 
        //MC_Step_Ising(int(randomNumber(0,Nx)),int(randomNumber(0,Ny)),T); 
	
	if(  (i%(stepsPerSite*Vol()/10))==0 ){ //every 10% output to terminal		 
		std::cout << ((double)i/Vol()) << " steps/site: E/vol = " 
		<< total_Energy() << ", " << E_total << ", P/vol = " 
                << total_Polarisation() << ", " << P_total << std::endl;
	}
        if(  (i%(stepsPerSite*Vol()/100))==0 ){ //every 1% output to file
            output <<(0.025*T)/(J*(double)300)<<" "<<(double)i/(Vol())<< " " 
            <<total_Energy()<<" "<< E_total<<" "<<total_Polarisation()<< "\n";	
        }
}
output.close();
}

//Lattice member func, runs statistics run on lattice
void Lattice::Run(int sampleDistance, int nSamples, double T) {
//Update global variables to latest values, they then have a starting 
//value for the run() function:
    Accepted=0; Rejected=0;
    E_total=total_Energy();
    Esqrd=E_total*E_total;
    P_total=total_Polarisation();
    Psqrd=P_total*P_total;
    
    //open output file;
    std::ofstream runOutput;
    runOutput.open("RunStats_"+std::to_string((int)T)+"K.dat");
    runOutput<<"#MCS E Esqrd E_av Esqrd_av P Psqrd P_sqrd_av\n";
    E_av=0.0;
    Esqrd_av=0.0;
    P_av=0.0;
    Psqrd_av=0.0;
    for (int j=0;j<=(nSamples);j++) {
    //We are in equilibrium so start running but updating the global
    //variables in the MC_step function. 
    for (int i=0;i<=(sampleDistance*Vol());i++) {
        //int acts like floor for numbers > 0 (but is ~x3 faster).
        MC_Step(int(randomNumber(0,Nx)),int(randomNumber(0,Ny)),T);
        //MC_Step_Ising(int(randomNumber(0,Nx)),int(randomNumber(0,Ny)),T);
    }
    //runOutput<<
    //Have been updating the estimators, now average them;
    //Not using continously updated estimators yet for checks ...
    E_av+=total_Energy();
    P_av+=total_Polarisation();
    Esqrd_av += total_Energy()*total_Energy(); //will use updated energy later, this is for checks first.
    Psqrd_av += total_Polarisation()*total_Polarisation();
    } //after nSamples take the average 
    //these are outputted in main.cpp loop.
    E_av=1.0*E_av/nSamples;
    P_av=1.0*P_av/nSamples;
    Esqrd_av=Esqrd_av/nSamples;
    Psqrd_av=Psqrd_av/nSamples;
    Cv=(Esqrd_av-E_av*E_av)/(Vol()*kB*T*T);
    Chi=Vol()*(Psqrd_av-P_av*P_av)/(kB*T);
}

void Lattice::MC_Step_Ising(int x, int y, double T){
    /* Calc sum of neighbouring spins (z-dir for now) */
    double sum=0.0;
    dipole p = get_dipole(x,y);
    double pz_current=p.z;
    //double beta=300/(0.025*T);
    
    for (int i=-1;i<=1;i+=2) {
	 p=get_dipole(x+i,y);
         sum += p.z;
	}
    for (int j=-1;j<=1;j+=2) {
	 p=get_dipole(x,y+j);
         sum += p.z;
	}
    /* Calc change in energy */
    double delta=2*J*sum*pz_current;
    if (delta<=0) {
        lattice[x+y*Nx].z = -lattice[x+y*Nx].z;
    }
    else if (randomNumber(0,1)<exp(-1.0*delta/T)) {
        lattice[x+y*Nx].z = -lattice[x+y*Nx].z;
    }

}

//Lattice member func, performs MC step on x,y,z coordinate dipole
void Lattice::MC_Step(int x, int y, double T) { 
/*Recover energy variable, should make that private?
Then propose random new dipole, calc change in energy for that new config
relative to old config, accept or reject new dipole.*/
//Note: see paper http://csml.northwestern.edu/resources/Preprints/mclr.pdf

//Save current dipole direction
dipole p_old = lattice[x+y*Nx]; //not useful for Ising but will be later.
//generate new trial dipole direction.
dipole p_new = move(x,y);
//energy of change if we flip dipole at (x,y).
double dE=deltaE(x,y, p_new); 

//If energy change is positive need Boltzmann factor vs rand number
if (dE>0) {
	if (exp(-dE/T)<randomNumber(0,1)) { //i.e if boltz loses, reject change.
                Rejected++;
                lattice[x+y*Nx]=p_old;
                //lattice[x+y*Nx].z=-lattice[x+y*Nx].z; //reverting to old
                 //change already occurred in deltaE().
		}
	else { //if move was accepted anyway update variables and keep change.
            Accepted++;
            E_total+=dE/Vol();
            Esqrd += (dE*dE)/(Vol()*Vol()); 
            P_total += 2*lattice[x+y*Nx].z/Vol(); 
            Psqrd += (2*lattice[x+y*Nx].z*2*lattice[x+y*Nx].z)/(Vol()*Vol());
        }
}
else if(dE<=0) {
//do nothing to updated dipole (updated in deltaE), 
// i.e accept the new move performed in deltaE(...).
    Accepted++;
    E_total += dE/Vol();
    Esqrd += (dE*dE)/(Vol()*Vol()); 
    //below are wrong for non-Ising z-models.
    P_total += 2*p_new.z/Vol(); //dP = 2s_new (see MC books e.g Newman and Barkema)
    Psqrd += (2*p_new.z*2*p_new.z)/(Vol()*Vol());
} 
}//end of MC_Step.

//Change in energy upon flip,
//has been abstracted to help
//statistics gathering.
double Lattice::deltaE(int x, int y, dipole p_new) {
//Calc current energy and dipole
double E_beforeFlip=site_Hamiltonian(x,y);//Energy(x,y);
//update lattice point as new dipole and calc new energy
lattice[x+y*Nx].x=p_new.x;
lattice[x+y*Nx].y=p_new.y;
lattice[x+y*Nx].z=p_new.z;
//std::cout << "lattice (px,py,pz)= " << lattice[x+y*Nx+z*Nx*Ny].x << " " << lattice[x+y*Nx+z*Nx*Ny].y << " " << lattice[x+y*Nx+z*Nx*Ny].z << std::endl;
double E_afterFlip=site_Hamiltonian(x,y);//Energy(x,y);
//double E_afterTEST=total_Energy();
//std::cout<< "Ediff = " << E_afterFlip-E_beforeFlip << std::endl;
double dE=(E_afterFlip-E_beforeFlip);
//Testing for errors.
/*if(dE != (E_afterTEST-E_beforeTEST)) {
std::cout<< "\nPROBLEM!!! (x,y)="<<x<<", "<<y<<"\ndE=" << dE << "\nE_after-E_before=" 
        << E_afterTEST-E_beforeTEST << "\ndP="<< p_new.z-p_old.z
        << "\nP_new="<<p_new.x<<", "<<p_new.y<<", "<<p_new.z
        << "\nP_old="<<p_old.x<<", "<<p_old.y<<", "<<p_old.z
        << "\nlattice[x,y]="<<lattice[x+y*Nx].x<<", "<<lattice[x+y*Nx].y
        <<", "<<lattice[x+y*Nx].z<<std::endl;
std::cout<<"\n Surrounding lattice:\n"
        << "    " << lattice[x+(y+1)*Nx].z << "    \n"
        << lattice[x-1+y*Nx].z << "   "<< lattice[x+y*Nx].z << "   " << lattice[x+1+y*Nx].z <<"\n"
        << "    " << lattice[x+(y-1)*Nx].z << "    \n" << std::endl;
*          
}*/

return dE;
}

Lattice::dipole Lattice::move(int x,int y) { //move dipole at 
    dipole p;
    if(model.compare("ISING")==0){
        p.x=lattice[x+y*Nx].x; //Ising so only z component.
        p.y=lattice[x+y*Ny].y; //Ising so only z component.
        p.z= -lattice[x+y*Nx].z; //Ising so can only flip.
    }
    if(model.compare("DIPOLE-DIPOLE")==0){
        //first just move in any direction, pi flips etc allowed.
        //list of triplets of the directions. Perhaps define externally when
        //I allow more directions.
        int dirs[18] = {-1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,1, 0,0,-1};
        int dir=int(randomNumber(0,6));
        p.x=dirs[dir*3]; //6 directions, each with x coordinates.
        p.y=dirs[dir*3+1];  
        p.z=dirs[dir*3+2]; 
    }
    return p;
}
/* This function is key for the different models*/
double Lattice::site_Energy(int x, int y) {
    double E=site_Hamiltonian(x,y); //need this as well.
    for (int i=-1;i<=1;i+=2) {
	E += site_Hamiltonian(x+i,y);
	}
	for (int j=-1;j<=1;j+=2) {
	E += site_Hamiltonian(x,y+j);		
	}
    return E;    
}

//Site Hamiltonian is interaction between particle at x and y, NOT the
//total energy of the site!!!! That is in site_Energy()
double Lattice::site_Hamiltonian(int x, int y) {
//calculate the Hamiltonian for a site with coord's x,y,z
	//ISING
	//must calc with three seperate loops so that we don't get
	//next to nearest neighbours from (x+1,y+1,k+1)
    double H=0.0;
    dipole pi = get_dipole(x,y);
    dipole pj;
    
    if (model.compare("ISING")==0) {
        for (int i=-1;i<=1;i+=2) {
            H += -J*dot_dipole(pi,get_dipole(x+i,y));
	}
	for (int j=-1;j<=1;j+=2) {
            H += -J*dot_dipole(pi,get_dipole(x,y+j));		
	}
    } //ISING ENERGY 
    
    if (model.compare("DIPOLE-DIPOLE")==0) {
        double rij=0.0;
        for (int i=-r_cut;i<=r_cut;i++) {
            for(int j=-r_cut;j<=r_cut;j++) {
                if(i==0 && j==0) continue; //or we'll get Nan's.
                //went outside of interaction range.
                if((i*i+j*j)>r_cut*r_cut) continue; 
                pj=get_dipole(x+i,y+j);
                rij=sqrt(i*i+j*j);
                //std::cout<<rij<<std::endl;
                H += dot_dipole(pi,pj)/(rij*rij*rij);
                //dot product manually to save casting (i,j) as vector
                H -= (pi.x*i+pi.y*j)*(pj.x*i+pj.y*j)/(rij*rij*rij);
            }
        }      
    } //DIPOLE-DIPOLE ENERGY
    return H;
}
double Lattice::total_Energy(){
double E=0.0;
 for(int i=0;i<Nx;i++) {
          for(int j=0;j<Ny;j++) {
                E+=site_Hamiltonian(i,j); 
          }
        }      
return E/Vol();
}
//Calc. total polarisation of lattice.
//Choice of words not sacrosanct, this also
//refers to the magnetisation for spin simulations.
double Lattice::total_Polarisation(){
//vector to store net pol.
dipole P_net;P_net.x=0.0;P_net.y=0.0;P_net.z=0.0;
	for(int i=0;i<Nx;i++) {
    	      for(int j=0;j<Ny;j++) {
					P_net.x += lattice[i+j*Nx].x;
					P_net.y += lattice[i+j*Nx].y;
					P_net.z += lattice[i+j*Nx].z;	
			}
		}
double P_net_mag=sqrt(dot_dipole(P_net,P_net));
return P_net_mag/Vol();
}
double Lattice::dot_dipole(Lattice::dipole p1, Lattice::dipole p2){
double dot=p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
return dot;
}
//create random number between min and max using mersenneTwister.
double Lattice::randomNumber(double min, double max){
        std::uniform_real_distribution<double> uniformDistribution(min, max);
        //std::cout << uniformDistribution(m_rng) << std::endl;
        return uniformDistribution(m_rng);
    }





















