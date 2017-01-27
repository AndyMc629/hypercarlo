/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "inputParser.h"

/*struct Config {
    int num;
    std::string str;
    double flt;
};*/
void loadConfig(Config& config) {
    std::ifstream fin("input.cfg");
    std::string line;
    while (getline(fin, line)) {
        std::istringstream sin(line.substr(line.find("=") + 1));
        if (line.find("nx") != -1)
            sin >> config.nx;
        else if(line.find("ny") !=-1)
            sin >> config.ny;
        else if (line.find("model") != -1)
            sin >> config.model;
        else if (line.find("T_min") != -1)
            sin >> config.T_min;
        else if (line.find("T_max") != -1)
            sin >> config.T_max;
        else if (line.find("dT") != -1)
            sin >> config.dT;
        else if (line.find("equilSteps") != -1)
            sin >> config.equilSteps;
        else if (line.find("ensembleSize") != -1)
            sin >> config.ensembleSize;
        else if (line.find("sampleFreq") != -1)
            sin >> config.sampleFreq;
        else if (line.find("initialLattice") != -1)
            sin >> config.initialLattice;
        else if (line.find("degOfFreedom") != -1)
            sin >> config.degOfFreedom;
        else if (line.find("rCutOff") != -1)
            sin >> config.rCutOff;
        else if (line.find("runType") != -1)
            sin >> config.runType;
        
    }
}