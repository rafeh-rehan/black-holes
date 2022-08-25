# black-holes
This repository consists of work done for my senior thesis in modified theories of gravity. 

Included is working code developed in FORTRAN 90 to compute properties and aspects of spherically symmetric black holes and results of the findings.  

## **black_hole_accretion**

This code is designed to calculate the percent differences of the peaks of flux, temperature, and luminosity distributions between two spherically symmetric geometries.

The Schwarzschild geometry is taken to be the reference case. Also provided are two conformally related geometries corresponding to scalar tensor theories of modified gravity. These conformally related geometries can be changed to study other spacetimes in comparison to the Schwarzschild solution.

Compile as: gfortran tools.f90 roots.f90 integrate.f90 structure.f90


## **black_hole_shadows**

This code is designed to draw shadows of spherically symmetric black holes. Given a metric and its derivatives, the geodesic equation can be solved symplectically in order to trace the path of light from a source location to an observer.

The Schwarzschild metric is used as the base case for spherically symmetric black holes and their asymptotic solutions. A second spherically symmetric metric found in https://arxiv.org/pdf/0806.0679.pdf is also included. This second metric can be changed to study other spacetimes in comparison with the Schwarzschild geometry.

Compile as: gfortran -O3 -fdefault-real-8 tools.f90 shadows.f90
