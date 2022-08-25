# black-hole-accretion
This code is designed to calculate the percent differences
of the peaks of flux, temperature, and luminosity distributions
between two spherically symmetric geometries. 
 
The Schwarzschild geometry is taken to be the reference case. 
Also provided are two conformally related geometries corresponding to 
scalar tensor theories of modified gravity. These conformally related 
geometries can be changed to study other spacetimes in comparison to the 
Schwarzschild solution. 
 
Compile as:
gfortran tools.f90 roots.f90 integrate.f90 structure.f90
