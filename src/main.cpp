/* Monte Carlo lattice model simulator, containing several model Hamiltonians.
 * 
 * As of 26/1/2017 these are:
 * 
 * 1) Basic Ising model.
 * 2) Dipolar model with n=6 orientations allowed.
 * 
 * in 2D and 3D. Main object is the Lattice class, which contains
 * all of the relevant methods needed for a Monte Carlo run (initialisation,
 * energy etc).
 * 
 * Input parameters are read from 'input.cfg' file.
 * Outputs summary statistics in 'Summary.dat' file.
 * 
 * Author: Andrew P. McMahon, Imperial College London.
 * Began: 25/1/2016

*/

#include<cstdlib>
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<random>
#include<ctime> //output times in log files
#include "Lattice.h"
#include "inputParser.h"

int main() {
    /********************************************************/ 
    /***************** READ IN RUN PARAMETERS ***************/ 
    /********************************************************/ 
    // current date/time based on current system
    time_t now = time(0);
    tm *ltm = localtime(&now);
    std::string datetime=std::to_string(ltm->tm_mday)+"_"+std::to_string(1+ltm->tm_mon)+"_"+std::to_string(1900+ltm->tm_year)+"_"+std::to_string(ltm->tm_hour)+std::to_string(ltm->tm_min);
    //std::cout<<"t = "<<datetime<<"\n";
    std::ofstream runSummary;
    runSummary.open("Summary_"+datetime+".dat");

    Config config;
    loadConfig(config);
    runSummary<<"#------------------------------------------------------\n"
        <<"#Config file loaded. Selected parameters:\n#\n" 
        <<"#Lattice = "<<config.nx<<"_"<<config.ny<<"\n"
        <<"#Model = "<<config.model<<" with n = "<<config.degOfFreedom<<"\n#\n"
        <<"#RunType = ANNEALING\n"
        <<"#r_cut = "<<config.rCutOff<<"\n"
        <<"#T_min, T_max, dT = "<<config.T_min<<", "<<config.T_max<<", "<<config.dT<<"\n#\n"
        <<"#equilSteps = "<<config.equilSteps<<"\n#ensembleSize = "<<config.ensembleSize
        <<"\n#------------------------------------------------------\n"
        << "#T(K) E_av Esqrd_av P_av Psqrd_av Cv Chi OP tau(px) tau(py) tau(pz) tau(E) tau(OP)\n"; 

    //INTIALISE LATTICE
    Lattice lattice=Lattice(config.nx,config.ny,config.model);
    lattice.degOfFreedom=config.degOfFreedom;
    lattice.r_cut=config.rCutOff;
    lattice.initialise_lattice(config.initialLattice);
    lattice.output_lattice("InitialState_"+datetime+".dat");

    
    /********************************************************/
    /********************** TEMP LOOP ***********************/
    /********************************************************/
    int T_counter=int((config.T_max-config.T_min)/config.dT);
    for (int i=0; i<=T_counter;i++) {
        //CURRENT TEMP
        double T=(i*config.dT)+config.T_min;
        //EQUILIBRATE
        lattice.Equilibrate(config.equilSteps,T);
        //OUTPUT EQUIL LATTICE
        std::string equilFile="Lattice_equil_"+std::to_string(config.equilSteps)+"_stepsPerSite_"+std::to_string(T)+"K.dat";
        lattice.output_lattice(equilFile);
        //PERFORM RUN
        lattice.Run(config.sampleFreq, config.ensembleSize,T);
   
        //OUTPUT TO TERMINAL
        std::cout << "T="<<T<<"K:\n"
        << "E_av="<<lattice.E_av << "\n"
        << "P_av="<<lattice.P_av << "\n"
        << "OP = "<<lattice.orderParam_av<<"\n"
        << "Accepted="<<lattice.Accepted<<" , Rejected="<<lattice.Rejected<< "\n\n";    
        //OUTPUT TO FILE
        runSummary << T << " " << lattice.E_av << " " << lattice.Esqrd_av 
        << " " <<lattice.P_av<< " " << lattice.Psqrd_av << " " << lattice.Cv 
        << " " << lattice.Chi<< " " <<lattice.orderParam_av<<" "<<lattice.tau_px<< " " <<lattice.tau_py
        << " " <<lattice.tau_pz<< " "<<lattice.tau_E<< " " <<lattice.tau_orderParam<<"\n";
    }//END OF TEMP LOOP
    //tidy up
    runSummary.close();
    return 0;
}
