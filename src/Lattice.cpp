/* Source code for the 'Lattice' class, which contains main bulk of run 
 * functions.
 */

#include "Lattice.h"
#include<fstream>
#include<random> //for MT algorithm.
#include "Constants.h" //constants like pi,kB etc.
#include<cmath>
#include "mathsFuncs.h" //for correct modulo function etc...
#include<list> // for autocorrelation lists.

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
        else if(s.compare("COL_ANTI_FERRO")==0) {
          std::cout << "COL ANTI FERRO CHOSEN" << std::endl;
          for(int i=0;i<Nx;i++) {
                for(int j=0;j<Ny;j++) {
                  lattice[i+j*Nx].x = 0.0;//std::pow(-1,i); //causes nice striping.
                  lattice[i+j*Nx].y = std::pow(-1,i);//0.0;
                  lattice[i+j*Nx].z = 0.0; //ensure norm=1.
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
                  lattice[i+j*Nx].z = sqrt(1-(lattice[i+j*Nx].x*lattice[i+j*Nx].x
                          +lattice[i+j*Nx].y*lattice[i+j*Nx].y) ); //ensure norm=1.
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
//will output data like x,y,z,px,py,px,E_site,V_site 
    std::ofstream output;
    output.open(datafile.c_str());
    for(int i=0;i<Nx;i++) {
        for(int j=0;j<Ny;j++) {
                output << i << " " << j << " " << "0 " << get_dipole(i,j).x << " "
                << get_dipole(i,j).y << " "
                << get_dipole(i,j).z << " " << site_Hamiltonian(i,j) 
                << " " << site_Potential(i,j) << "\n";
                }
            }
    output.close();
}

//Lattice member func, equilibrate for stepsPerSite MCSweeps.
void Lattice::Equilibrate(int stepsPerSite, double T) {
//initialise the global var's.
E_total=total_Energy();
P_total=total_Polarisation();
P_total.norm=sqrt(dot_dipole(P_total,P_total));

std::cout << "Initial lattice E = " << E_total << "\n"
          << "Initial lattice P = (" 
          << P_total.x <<","<<P_total.y<<","<<P_total.z<<")\n";

//Open data file to store equilibration stats.
std::ofstream output;
output.open("EquilibrationStats_"+std::to_string(T)+"K.dat");
output << "#Equilibrated stats.\n"
        << "#kT/J MCS E_tot/site E_tot(global)/site P_tot/site\n";   

for (int i=0;i<=(stepsPerSite*Vol());i++) {
    if(  (i%(stepsPerSite*Vol()/10))==0 ){ //every 10% output to terminal
            dipole pp=total_Polarisation();
		std::cout << ((double)i/Vol()) << " steps/site: E/vol = " 
		<< total_Energy() << ", " << E_total << ", P/vol = " 
                << sqrt(dot_dipole(pp,pp)) << ", " << P_total.norm << std::endl;
	}
        if(  (i%(stepsPerSite*Vol()/100))==0 ){ //every 1% output to file
            dipole pp=total_Polarisation();
            output <<(0.025*T)/(J*(double)300)<<" "<<(double)i/(Vol())<< " " 
            <<total_Energy()<<" "<< E_total<<" "<<sqrt(dot_dipole(pp,pp))<< "\n";	
        }
        MC_Step(int(randomNumber(0,Nx)),int(randomNumber(0,Ny)),T); 
        //MC_Step_Ising(int(randomNumber(0,Nx)),int(randomNumber(0,Ny)),T); 		
}
output.close();
}

//Lattice member func, runs statistics run on lattice
void Lattice::Run(int sampleFreq, int nSamples, double T) {
    /*************************************************************/
    /************************ INITIALISE *************************/
    /*************************************************************/
    /*Update global variables to latest values, they then have a starting 
     * value for the run() function:*/
    Accepted=0; Rejected=0;
    E_total=total_Energy();
    Esqrd=E_total*E_total;
    P_total=total_Polarisation();
    Psqrd=P_total.norm*P_total.norm;
    orderParam_total=OrderParameter();
    /**************** AVERAGE VARIABLES (GLOBAL) *****************/
    E_av=0.0; //total 
    Esqrd_av=0.0;
    P_av=0.0;
    Psqrd_av=0.0;
    orderParam_av=0.0;
    int i,j; 
    /*************************************************************/
    /*************** INITIALISE CORRELATION **********************/
    /*************************************************************/
    
    /************ AUTOCORRELATION VARIABLES (GLOBAL) *************/
    PxCorr_av=0.0, PyCorr_av=0.0, PzCorr_av=0; //for autocorrelation Ising tests.
    ECorr_av=0.0;
    orderParamCorr_av=0.0;
    /************ AUTOCORRELATION VARIABLES (LOCAL) *************/
    //http://www.physics.buffalo.edu/phy410-505/2011/topic5/app2/
    const int nSave=100; //how many values to save for autocorr, make this 10 and sim dies ?
    std::list<double> pxSave,pySave,pzSave,ESave,OPSave;    
    //polarisation/magnetisation accumulators
    double px_autocorr[nSave]={0}, py_autocorr[nSave]={0}, pz_autocorr[nSave]={0};
    double E_autocorr[nSave]={0}, orderParam_autocorr[nSave]={0};
    int nCorr=0;
    //open output file;
    std::ofstream runOutput;
    runOutput.open("RunStats_"+std::to_string(T)+"K.dat");
    runOutput<<"# Autocorrelation functions of each variable listed.\n"
             <<"# Values for orderParam only listed if non-Ising model.\n"
             <<"# MCS Px Py Pz E OrderParam\n";
    /*************************************************************/
    /******************* PERFORM MC SWEEPS ***********************/
    /*************************************************************/
    for (j=0;j<=(nSamples);j++) {
    //We are in equilibrium so start running but updating the global
    //variables in the MC_step function. 
    for (i=0;i<=(sampleFreq*Vol());i++) {
        //int acts like floor for numbers > 0 (but is ~x3 faster).
        MC_Step(int(randomNumber(0,Nx)),int(randomNumber(0,Ny)),T);
        //MC_Step_Ising(int(randomNumber(0,Nx)),int(randomNumber(0,Ny)),T);        
    }   
    /*************************************************************/
    /******************* UPDATE ESTIMATORS  **********************/
    /*************************************************************/
    //Have been updating the estimators, now average them;
    //Not using continously updated estimators yet for checks ...
    E_av+=E_total; // safe to do this now if needs checking return to total_Energy();
    P_av+=sqrt(dot_dipole(total_Polarisation(),total_Polarisation())); //two calls = horrendous!
    Esqrd_av += total_Energy()*total_Energy(); //will use updated energy later, this is for checks first.
    Psqrd_av += dot_dipole(total_Polarisation(),total_Polarisation());
    orderParam_total = OrderParameter();
    orderParam_av += orderParam_total;
    
    /*************************************************************/
    /***********UPDATE AUTOCORRELATION ACCUMULATORS***************/
    /*************************************************************/
    //Update autocorrelation array
    if(pxSave.size()==nSave) {
        nCorr++;
        px_autocorr[0]+=P_total.x*P_total.x;
        py_autocorr[0]+=P_total.y*P_total.y;
        pz_autocorr[0]+=P_total.z*P_total.z;
        E_autocorr[0]+=E_total*E_total;
        orderParam_autocorr[0]+=orderParam_total*orderParam_total;
        PxCorr_av+=P_total.x;PyCorr_av+=P_total.y;PzCorr_av+=P_total.z;
        ECorr_av+=E_total;orderParamCorr_av+=orderParam_total;
        
        std::list<double>::const_iterator ipx = pxSave.begin(),ipy = pySave.begin(),ipz=pzSave.begin();
        std::list<double>::const_iterator ipE = ESave.begin(), ipOP = OPSave.begin();
    for (int ii = 1; ii <= nSave; ii++) {
        px_autocorr[ii] += *ipx++ * P_total.x;
        py_autocorr[ii] += *ipy++ * P_total.y;
        pz_autocorr[ii] += *ipz++ * P_total.z;
        E_autocorr[ii] += *ipE++ * E_total;
        orderParam_autocorr[ii] += *ipOP++ * orderParam_total;
    }
        //discard oldest values
        pxSave.pop_back();pySave.pop_back();pzSave.pop_back();
        ESave.pop_back();OPSave.pop_back();
    }
    //save current values
    pxSave.push_front(P_total.x);pySave.push_front(P_total.y);pzSave.push_front(P_total.z);
    ESave.push_front(E_total);OPSave.push_front(orderParam_total);
    } //after nSamples take the average 
    /*************************************************************/
    /******************* AVERAGE ESTIMATORS **********************/
    /*************************************************************/
    //these are outputted in main.cpp loop.
    E_av=1.0*E_av/nSamples;
    P_av=1.0*P_av/nSamples;
    Esqrd_av=Esqrd_av/nSamples;
    Psqrd_av=Psqrd_av/nSamples;
    Cv=(Esqrd_av-E_av*E_av)/(Vol()*kB*T*T);
    Chi=Vol()*(Psqrd_av-P_av*P_av)/(kB*T);
    orderParam_av=orderParam_av/nSamples;
    /*************************************************************/
    /*********** OUTPUT AUTOCORRELATION FUNCS ********************/
    /*************************************************************/    
    //calculate autocorrelation time for Pz, Pz_AutoCorrTime is a 
    //public var in the 'Lattice' class
    //double avpxCorr=PxCorr_av/nCorr,avpyCorr=PyCorr_av/nCorr,avpzCorr=PzCorr_av/nCorr;
    //double avECorr=ECorr_av/nCorr, avOPCorr=orderParamCorr_av/nCorr;
    PxCorr_av=PxCorr_av/nCorr;PyCorr_av=PyCorr_av/nCorr;PxCorr_av=PzCorr_av/nCorr;
    ECorr_av=ECorr_av/nCorr; orderParamCorr_av=orderParamCorr_av/nCorr;
    
    double cpx0 = px_autocorr[0]/nCorr - PxCorr_av*PxCorr_av;
    double cpy0 = py_autocorr[0]/nCorr - PyCorr_av*PyCorr_av;
    double cpz0 = pz_autocorr[0]/nCorr - PzCorr_av*PzCorr_av;
    double cE0 = E_autocorr[0]/nCorr - ECorr_av*ECorr_av;
    double cOP0 = orderParam_autocorr[0]/nCorr - orderParamCorr_av*orderParamCorr_av;
    //re-initialise the autocorrelation times
    tau_px=0.0, tau_py=0.0, tau_pz=0, tau_E=0.0, tau_orderParam=0.0;
    for (int k=1;k<=nSave;k++) {
        runOutput<<k<<" "<<(px_autocorr[k]/nCorr-PxCorr_av*PxCorr_av)/cpx0
                <<" "<<(py_autocorr[k]/nCorr-PyCorr_av*PyCorr_av)/cpy0
                <<" "<<(pz_autocorr[k]/nCorr-PzCorr_av*PzCorr_av)/cpz0
                <<" "<<(E_autocorr[k]/nCorr-ECorr_av*ECorr_av)/cE0
                <<" "<<(orderParam_autocorr[k]/nCorr-orderParamCorr_av*orderParamCorr_av)/cOP0<<"\n";
        tau_px += (px_autocorr[k]/nCorr-PxCorr_av*PxCorr_av)/cpx0;
        tau_py += (py_autocorr[k]/nCorr-PyCorr_av*PyCorr_av)/cpy0;
        tau_pz += (pz_autocorr[k]/nCorr-PzCorr_av*PzCorr_av)/cpz0;
        tau_E += (E_autocorr[k]/nCorr-ECorr_av*ECorr_av)/cE0;
        tau_orderParam += (orderParam_autocorr[k]/nCorr-orderParamCorr_av*orderParamCorr_av)/cOP0;
    }  
    runOutput.close();
}

double Lattice::OrderParameter() {
    double OP=0.0; //order param
    if(model.compare("DIPOLE-DIPOLE")==0){
        /* See Tomita, Monte Carlo Study of Two-Dmensional Heisenberg
         * Dipolar Lattices.
         * Journal of the Physical Society of Japan Vol. 78, 
         * No. 11, November, 2009, 114004 
         */
        // for dipole-dipole model OP is a vector 
        double OPx=0.0,OPy=0.0,OPz=0.0;
        dipole p;
        for(int i=0;i<Nx;i++) {
            for(int j=0;j<Ny;j++) {
                p=get_dipole(i,j);
                OPx+=std::pow(-1,j)*p.x;//std::pow(-1,mod(j,2))*p.x;
                //std::cout<<"i,j ="<<i<<", "<<j<<" px="<<p.x<<"\n";
                OPy+=std::pow(-1,i)*p.y;//std::pow(-1,mod(i,2))*p.y;
                OPx+=std::pow(-1,i+j)*p.z;//std::pow(-1,mod(i+j,2))*p.z;
            }
        }
        OP=OPx*OPx+OPy*OPy+OPz*OPz;
        //std::cout<< "OP = "<<OP<<"\n";
    }
    return OP;
}

void Lattice::MC_Step_Ising(int x, int y, double T){
    /* Calc sum of neighbouring spins (z-dir for now) */
    double sum=0.0;
    dipole p = get_dipole(x,y);
    double pz_current=p.z;
    //std::cout << "E_total in MC(before)= "<<E_total<<", "<<total_Energy()<<"\n";
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
        E_total+=1.0*delta;
    }
    else if (randomNumber(0,1)<exp(-1.0*delta/T)) {
        lattice[x+y*Nx].z = -lattice[x+y*Nx].z;
        E_total+=1.0*delta;
    }
    //std::cout<<"delta = "<<delta<<"\n"
    //<< "E_total in MC(after)= "<<E_total<<", "<<total_Energy()<<"\n\n";
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
            E_total+=dE;///Vol();
            Esqrd += (dE*dE);///(Vol()*Vol()); 
            P_total.x+=(p_new.x-p_old.x);//lattice[x+y*Nx].x;///Vol();
            P_total.y+=(p_new.y-p_old.y);//lattice[x+y*Nx].y;///Vol();
            P_total.z+=(p_new.z-p_old.z);//lattice[x+y*Nx].z;///Vol();
            P_total.norm=sqrt(dot_dipole(P_total,P_total));
            Psqrd += P_total.norm*P_total.norm;
        }
}
else if(dE<=0) {
//do nothing to updated dipole (updated in deltaE), 
// i.e accept the new move performed in deltaE(...).
    Accepted++;
    E_total += dE;//Vol();
    Esqrd += (dE*dE);///(Vol()*Vol()); //will divide by vol at end.
    //below are wrong for non-Ising z-models.
    P_total.x+=(p_new.x-p_old.x);//lattice[x+y*Nx].x;///Vol(); // INCORRECT, NEED TO ADD THE CHANGE NOT THE VALUE!
    P_total.y+=(p_new.y-p_old.y);//lattice[x+y*Nx].y;///Vol();
    P_total.z+=(p_new.z-p_old.z);//lattice[x+y*Nx].z;///Vol();
    P_total.norm=sqrt(dot_dipole(P_total,P_total));
    Psqrd += P_total.norm*P_total.norm;
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
    if(model.compare("DIPOLE-DIPOLE")==0){// && degOfFreedom==6){
        //first just move in any direction, pi flips etc allowed.
        //list of triplets of the directions. Perhaps define externally when
        //I allow more directions.
        //int dirs[18] = {-1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,1, 0,0,-1};
    //theta, phi possibilities
    double dirs[12]={0,0, PI,0, 
                    0.5*PI,0, 0.5*PI,0.5*PI, 0.5*PI,PI, 0.5*PI,1.5*PI};
    
        int dir=int(randomNumber(0,6)); //=faces of cube.
        double theta=dirs[dir*2];
        double phi=dirs[dir*2+1];
        p.x=sin(theta)*cos(phi);//dirs[dir*3]; //6 directions, each with x coordinates.
        p.y=sin(theta)*sin(phi);//dirs[dir*3+1];  
        p.z=cos(theta);//dirs[dir*3+2]; 
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
                H -= (pi.x*i+pi.y*j)*(pj.x*i+pj.y*j)/(rij*rij*rij*rij*rij);
            }
        }      
    } //DIPOLE-DIPOLE ENERGY
    return H;
}

double Lattice::site_Potential(int x, int y){
    double Vxy=0.0;
    dipole pj;
    if (model.compare("DIPOLE-DIPOLE")==0) {
        double rij=0.0;
        for (int i=-r_cut;i<=r_cut;i++) {
            for(int j=-r_cut;j<=r_cut;j++) {
                if(i==0 && j==0) continue; //or we'll get Nan's.
                //went outside of interaction range.
                if((i*i+j*j)>r_cut*r_cut) continue; 
                pj=get_dipole(x+i,y+j);
                rij=sqrt(i*i+j*j);
                //dot product manually to save casting (i,j) as vector
                Vxy += (pj.x*i+pj.y*j)/(rij*rij);
            }
        }      
    } //DIPOLE-DIPOLE POTENTIAL
    return Vxy;
}
double Lattice::total_Energy(){
double E=0.0;
 for(int i=0;i<Nx;i++) {
          for(int j=0;j<Ny;j++) {
                E+=0.5*site_Hamiltonian(i,j); //0.5 to stop double counting??? Makes it work.
          }
        }      
return E;///Vol();
}
//Calc. total polarisation of lattice.
//Choice of words not sacrosanct, this also
//refers to the magnetisation for spin simulations.
Lattice::dipole Lattice::total_Polarisation(){
//vector to store net pol.
dipole P_net;P_net.x=0.0;P_net.y=0.0;P_net.z=0.0;
	for(int i=0;i<Nx;i++) {
    	      for(int j=0;j<Ny;j++) {
					P_net.x += lattice[i+j*Nx].x;
					P_net.y += lattice[i+j*Nx].y;
					P_net.z += lattice[i+j*Nx].z;	
			}
		}
//double P_net_mag=sqrt(dot_dipole(P_net,P_net));
//will divide by Vol() when outputting? same time as divisions take more time?
//P_net.x*=P_net.x/Vol(); P_net.y*=P_net.y/Vol(); P_net.z*=P_net.z/Vol();
return P_net;
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




















