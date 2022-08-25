
This code is designed to draw shadows of spherically symmetric 
black holes. Given a metric and its derivatives, the geodesic equation 
can be solved symplectically in order to trace the path of light from 
a source location to an observer. 

The Schwarzschild metric is used as the base case for spherically 
symmetric black holes and their asymptotic solutions. A second 
spherically symmetric metric found in https://arxiv.org/pdf/0806.0679.pdf
is also included. This second metric can be changed to study other 
spacetimes in comparison with the Schwarzschild geometry.

An implicit Runge-Kutta method for solving sets of differential equations
is used to solve the geodesic equations (reference: https://arxiv.org/abs/1704.04114)
 
Compile as:
gfortran -O3 -fdefault-real-8 tools.f90 shadows.f90

 
