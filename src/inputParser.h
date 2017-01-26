/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   inputParser.h
 * Author: andrew
 *
 * Created on 26 January 2017, 12:24
 */

#ifndef INPUTPARSER_H
#define INPUTPARSER_H
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>

struct Config{
    int nx; int ny; 
    std::string model; 
    int degOfFreedom; //only applies to dipolar models for now
    int rCutOff;
    double T_min; double T_max; double dT;
    int equilSteps; int ensembleSize;
    int sampleFreq;
    std::string initialLattice;
};
void loadConfig(Config&);


#endif /* INPUTPARSER_H */

