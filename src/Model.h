/* 
 * File:   Model.h
 * Author: apm13
 *
 * Created on 17 January 2017, 13:18
 * 
 * Header file for the "Model" class.
 */

#ifndef MODEL_H
#define MODEL_H
#include<string>
#include<iostream>
#include "Lattice.h"

class Lattice; //Model is friend of Lattice so need forward declaration.
class Model {
    friend class Lattice;
public:
    //Model constructor
    Model(std::string);
    //Model default constructor
    Model();
    //Model destructor
    ~Model(); 
    double Ising_Site_Hamiltonian(int, int);
private:
    std::string modelName;  
protected:
    Lattice lattice;
};
#endif /* MODEL_H */

