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
#include<ctime> //output times in log files
//using namespace std;

int main() {

 /*************************** 
 * SET UP THE RUN PARAMETERS 
 * will do this from a config file eventually. 
 ****************************/    
//Set size of lattice
int Nx=10,Ny=10;//Nx=40,Ny=40;
//Seed a Mersenne Twister random number generator with current time.
std::mt19937 rng{static_cast<std::mt19937>(std::chrono::high_resolution_clock::now().time_since_epoch().count())};

//Initialise a lattice object.
Lattice lattice=Lattice(Nx,Ny,rng);

//Initialise the lattice to FERRO or PARA_ISING for Ising model.
lattice.initialise_lattice("FERRO");
//output initial lattice.
lattice.output_lattice("InitialState.dat");

//equilibrate lattice at temp T.
//double temp; //K
double T_min=0.2;//1000; //K
double T_max=5.0;//2000; //K
double dT=0.2; //K
int T_counter=int((T_max-T_min)/dT);
int equilStepsPerSite=2000;//10000;//100000;//10000;
int ensemble_size=1000000;

std::ofstream mainOutput;
mainOutput.open("Output.dat");

// current date/time based on current system
time_t now = time(0);
// convert now to string form
char* dt = ctime(&now);

mainOutput << "# Run began:" << dt << "\n" 
           << "# Key run parameters:\n"
           << "# Nx Ny = "<<Nx << " "<<Ny<< "\n"
           << "# (T_min,T_max,dT) = ("<<T_min<<", "<<T_max<<", "<<dT<<")\n"
           << "# equilStepsPerSite = " << equilStepsPerSite << "\n"
           << "# ensemble_size = " << ensemble_size << "\n#\n";
mainOutput << "#T(K) E_av Esqrd_av P_av Psqrd_av Cv Chi\n"; 

//for (int T=T_min;T<=T_max;T=T+dT) {
for (int i=0; i<=T_counter;i++) {
double T=(i*dT)+T_min;

//If you want to re-initialise at each temp (bad idea generally).
//lattice.initialise_lattice("FERRO");
//equilibrate
lattice.Equilibrate(equilStepsPerSite,T);
//output equilibrated lattice.
std::string equilFile="Lattice_equil_" + std::to_string(equilStepsPerSite)+"_stepsPerSite_"+std::to_string(T)+"K.dat";
lattice.output_lattice(equilFile);
//calc estimators
lattice.Run(1, ensemble_size,T);

std::cout << "T="<<T<<"K:\n"
<< "E_av="<<lattice.E_av << "\n"
<< "P_av="<<lattice.P_av << "\n"
<< "Accepted="<<lattice.Accepted<<" , Rejected="<<lattice.Rejected<< "\n\n";

mainOutput << T << " " << lattice.E_av << " " << lattice.Esqrd_av << " "
<<lattice.P_av<< " " << lattice.Psqrd_av << " " << lattice.Cv 
<< " " << lattice.Chi<<"\n";

}
mainOutput.close();
return 0;
}
