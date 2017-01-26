H  H Y   Y PPPP EEEE RRRR CCCC AAAAA RRRR L     OOOOO  
H  H  Y Y  P  P E    R  R C    A   A R  R L     O   O  
HHHH   Y   PPPP EEEE RRRR C    AAAAA RRRR L     O   O  
H  H   Y   P    E    R R  C    A   A R R  L     O   O  
H  H   Y   P    EEEE R  R CCCC A   A R  R LLLLL OOOOO    

# HyPerCarlo 

HyPerCarlo is a Monte Carlo lattice simulator with some models of use
for researching Hybrid Perovskites. 

Such models include (in 2D and 3D versions of the code):

1) Ising model  
2) Various n-Potts models with Ising like interactions  
3) Random field Ising model  
4) Random field n-Potts models  
5) Dipolar lattice model (with n-Potts constrained orientations)  
6) Dipolar lattice model (with random fields).  

# Background

The initial work on this code was inspired by Jarvist Frost's StarryNight code, published in:

"Atomistic Origins of High-Performance in Hybrid Halide Perovskite Solar Cells".
Nano Lett., 2014, 14 (5), pp 2584–2590, http://pubs.acs.org/doi/abs/10.1021/nl500390f.

and used by myself in:

"The dynamics of methylammonium cations in hybrid organic-inorganic perovskite solar cells"
Nature Communications 6, Article no 7124, (20150), http://www.nature.com/articles/ncomms8124
 
However, the current code is slightly more flexible through its inclusion of several more models,
which I am performing simulations on for my PhD Thesis.

This code will also use an Ewald summation for the calculation of long-range forces as an option.

I also want to eventually let users specify their own lattice Hamiltonians in input files, however 
in initial stages these will be hard-coded in. 

Note to self: Possible Ewald algorithm
http://csml.northwestern.edu/resources/Preprints/mclr.pdf.

Andrew P. McMahon, 
Theory and Simulation of Materials, 
Imperial College London.

Code began 25/1/2016.
