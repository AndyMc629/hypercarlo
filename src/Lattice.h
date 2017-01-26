#ifndef LATTICE_H
#define LATTICE_H
#include<vector>
#include<string>
#include<iostream>
#include<random>
#include<chrono>
#include<list> // for autocorrelation lists.

class Lattice {
    //friend class Model;
public:
    /*************************************************************/
    /************************* FUNCTIONS *************************/
    /*************************************************************/
    //Initialise a lattice/lattice constructor 
    //Lattice (int,int,std::mt19937&);
    Lattice(int,int, std::string);
    Lattice();
    //Lattice destructor
    ~Lattice();
    //Grabs volume of lattice/crystal.
    int Vol(void);
    /*Define dipole struct. NB: next 2 lines 
     * were in private before but I want user to
     * be able to access these, I think.
     * i.e means I can change some vectors ad hoc.*/
    struct dipole {double x; double y; double z; double norm;};
    //the lattice object is a vector of dipole structs.
    std::vector<dipole> lattice;	
    //get the dipole at xyz.
    dipole get_dipole(int,int);
    dipole move(int,int);
    /*initialisation state dependent on string
     *  string can have values:
     * 	FERRO = ferro arrangement 1,0,0 or something
     *  PARA = paraelectric , rand() ...
     * 	PREV = read in previous state as initial state
     * 	NB: I'll prob make this private later ...?*/
    void initialise_lattice(std::string);
    /*Output to file supplied by string*/
    void output_lattice(std::string);
    void Equilibrate(int, double);
    void Run(int,int,double);
    double OrderParameter();
    double site_Hamiltonian(int, int);
    double site_Potential(int, int);
    double site_Energy(int, int);
    double dot_dipole(Lattice::dipole p1, Lattice::dipole p2);
    double randomNumber(double, double);
    //std::mt19937& m_rng;
    //void initialise_lattice(std::string);
    void MC_Step(int,int,double);
    void MC_Step_Ising(int,int,double);
    double total_Energy();
    dipole total_Polarisation();
    double deltaE(int,int,Lattice::dipole);
    /*************************************************************/
    /************************* VARIABLES *************************/
    /*************************************************************/
    double E_total=0.0,Esqrd=0.0,Psqrd=0.0;
    dipole P_total;// P_total.x=0.0,P_total.y=0.0,P_total.z=0.0; //need this for correct calc.
    double E_av=0.0,P_av=0.0,Esqrd_av=0.0,Psqrd_av=0.0,Cv=0.0,Chi=0.0;
    double orderParam_total=0.0, orderParam_av=0.0; //will have a diff order param for different models
    /* Autcorrelation variables*/
    double ECorr_av=0.0;
    double PxCorr_av=0.0, PyCorr_av=0.0, PzCorr_av=0.0; 
    double orderParamCorr_av=0.0;
    double tau_E, tau_px, tau_py, tau_pz, tau_orderParam; //autocorrelation times
    //For MC stats.
    int Accepted=0,Rejected=0;
    /* Length of interactions in Dipole-Dipole models.
    * Public so's I can output it in main during testing but 
    * would rather it private. Will make a function called "run_summarise()" 
    * that outputs the relevant info. */
    int r_cut=2; 
    int degOfFreedom; //degrees of freedom for dipole models
    /*************************************************************/
private:        
    //sizes of crystal
    int Nx,Ny;
    double J=1.0;//J=0.025;
    //int r_cut=2;//1; //cut off for dipole-dipole interaction
    //seed Mersenne Twister algorithm.
    std::mt19937 m_rng{static_cast<std::mt19937>(std::chrono::high_resolution_clock::now().time_since_epoch().count())};
    std::string model;
};

#endif
