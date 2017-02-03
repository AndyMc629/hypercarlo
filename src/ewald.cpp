/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.

 * 
 * 
 * SEE "Computer simulations of dipolar fluids using ewald summations" in my 
 * Ewald SUmmation Mendeley folder and compare with ewald.cpp in "EwaldSummation"
 * folder on GitHub (pretty certain they are the same so can use this! I think...,
 * be careful to make sure that one is for dipoles and not for charges!
 *  */
#include<cmath>
#include "ewald.h"
#include "Lattice.h"
#include "Constants.h"
using namespace Constants;
// Computes the real and recriprocal summation in Ewald summation
double RealandReciprocalSpace(Lattice &lattice) 
{
    std::cout <<"Real and reciprocal is running!\n";
    
    //Initial values
    double V = 1.0;//(vol of one cell?)
    double Lx=1.0, Ly=1.0,Lz=1.0; //will allow me to change unit cell dims.
    double u1=0.0,u2=0.0;
    double xx = 0, yy = 0, zz = 0; //rij components.
    int acut=lattice.r_cut, gcut=5; //cut off in realc and reciprocal space.
    
    //Reciprocal vectors
    //Reciprocal space variable
    double g1, g2, g3, g_2, g_x, g_y, g_z;
    g1 = (2 * PI / V) * (Ly * Lz);
    g2 = (2 * PI / V) * (Lz * Lx);
    g3 = (2 * PI / V) * (Lx * Ly);
        
    return 0; 
}
    
    /*
 double V = Lx * Ly * Lz;
 double u1 = 0.0, u2 = 0.0;
 double xx = 0, yy = 0, zz = 0;
 int acut = 50.0, gcut = 50.0;

  //Reciprocal space variable
  double g1, g2, g3, g_2, g_x, g_y, g_z;
  g1 = (2 * Pi / V) * (Ly * Lz);
  g2 = (2 * Pi / V) * (Lz * Lx);
  g3 = (2 * Pi / V) * (Lx * Ly);
  
  //Real space variable
  double a1 = 0, a2 = 0, a3 = 0, len = 0;
  for(unsigned int i = 0;i < r.size();i += 4)
    {
      for(unsigned int j = 0;j < r.size();j += 4)
        {
          xx = r[j+0] - r[i+0];
          yy = r[j+1] - r[i+1];
          zz = r[j+2] - r[i+2];
      
          for(int x = -ds;x <= ds;x++)
            {
              for(int y = -ds;y <= ds;y++)
                {
                  for(int z = -ds;z <= ds;z++)
                    {  
                       //Reciprocal Space
                       g_x = g1*x;
                       g_y = g2*y;
                       g_z = g3*z;
                       if(norm(g_x,g_y,g_z) < gcut)
                       {
                         if(x==0&&y==0&&z==0){}
                         else
                         {
                           g_2 = dp(g_x,g_y,g_z,g_x,g_y,g_z);
                           u1 +=  (r[i+3] * r[j+3])*(1/g_2)*exp(-g_2/(4.0*kappa*kappa))*cos(dp(g_x,g_y,g_z,xx,yy,zz));
                         }
                       }

                       //Real Space
                       a1 = x * Lx + xx;
                       a2 = y * Ly + yy;
                       a3 = z * Lz + zz;
                       len = norm(a1,a2,a3); 
                       if(norm(x*Lx,y*Ly,z*Lz) < acut)
                       {
                         if( x==0 && y==0 && z==0)
                         {
                           if(i==j){}
                           else
                           { 
                             u2 +=  (r[i+3] * r[j+3]) * erfc(kappa*len)/len;
                           }
                         }
                         else
                         {
                           u2 += (r[i+3] * r[j+3]) * erfc(kappa*len)/len;
                         }
                         
                       }

                    }       
                }
            }
        }
    }

  return (u1*4.0*Pi * 0.5/V) + (u2/2.0);
} */
/*
// Computes the point energy
double PointEnergy(vector<double> const &r)
{
  double u = 0;
  for(unsigned int i = 0;i < r.size();i += 4)
    { 
      u = u - r[i+3]*r[i+3];
    }
  return u;
}

// Computes the Dipole
double Dipole(vector<double> const &r)
{
  double P = 0, P_x = 0,P_y = 0,P_z = 0;
  for(unsigned int i = 0; i < r.size();i += 4)
    {
      P_x += r[i+0]*r[i+3];
      P_y += r[i+1]*r[i+3];
      P_z += r[i+2]*r[i+3];
    }
  return  norm(P_x,P_y,P_z);
}*/